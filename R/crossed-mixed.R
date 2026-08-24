.canonical_term <- function(x) paste(sort(strsplit(x, ":", fixed = TRUE)[[1L]]), collapse = ":")

.pseudo_inverse_solve <- function(A, b, tol = 1e-10) {
  decomposition <- svd(A)
  keep <- decomposition$d > tol * max(decomposition$d)
  if (!any(keep)) return(rep(0, ncol(A)))
  as.vector(decomposition$v[, keep, drop = FALSE] %*%
    ((t(decomposition$u[, keep, drop = FALSE]) %*% b) / decomposition$d[keep]))
}

#' Construct an EMS-matched quasi-denominator
#'
#' Finds coefficients for candidate mean squares whose combined EMS matches
#' the numerator EMS. An exact single mean square is preferred; otherwise a
#' minimum-norm linear combination is returned with Satterthwaite degrees of
#' freedom.
#'
#' @param numerator_ems Named numeric EMS coefficient vector.
#' @param denominator_ems Matrix with one candidate mean square per row and EMS
#'   components in columns.
#' @param mean_squares Candidate observed mean squares.
#' @param dfs Candidate degrees of freedom.
#' @param tol Numerical matching tolerance.
#' @return A machine-readable denominator specification.
#' @export
quasi_denominator <- function(numerator_ems, denominator_ems, mean_squares, dfs,
                              tol = 1e-10) {
  denominator_ems <- as.matrix(denominator_ems)
  if (!is.numeric(numerator_ems) || !is.numeric(denominator_ems)) {
    stop("EMS coefficients must be numeric.", call. = FALSE)
  }
  if (ncol(denominator_ems) != length(numerator_ems)) {
    stop("Candidate and numerator EMS components do not conform.", call. = FALSE)
  }
  if (nrow(denominator_ems) != length(mean_squares) || length(dfs) != length(mean_squares)) {
    stop("One mean square and degrees of freedom value are required per candidate.", call. = FALSE)
  }
  if (any(!is.finite(mean_squares)) || any(!is.finite(dfs)) || any(dfs <= 0)) {
    stop("Candidate mean squares and degrees of freedom must be finite, with positive df.", call. = FALSE)
  }

  differences <- apply(denominator_ems, 1L, function(row) max(abs(row - numerator_ems)))
  exact_rows <- which(differences <= tol)
  if (length(exact_rows)) {
    coefficients <- rep(0, nrow(denominator_ems))
    coefficients[exact_rows] <- dfs[exact_rows] / sum(dfs[exact_rows])
    exact_match <- TRUE
    exact_single <- length(exact_rows) == 1L
  } else {
    coefficients <- .pseudo_inverse_solve(t(denominator_ems), numerator_ems, tol = tol)
    mismatch <- max(abs(as.vector(t(denominator_ems) %*% coefficients) - numerator_ems))
    if (!is.finite(mismatch) || mismatch > sqrt(tol)) {
      stop("No EMS-matched denominator can be constructed from the candidate mean squares.", call. = FALSE)
    }
    coefficients[abs(coefficients) <= tol] <- 0
    exact_match <- FALSE
    exact_single <- FALSE
  }

  value <- sum(coefficients * mean_squares)
  variance_term <- sum((coefficients * mean_squares)^2 / dfs)
  denominator_df <- if (exact_match) {
    sum(dfs[exact_rows])
  } else if (variance_term > tol) value^2 / variance_term else NA_real_
  names(coefficients) <- rownames(denominator_ems)

  list(
    coefficients = coefficients,
    value = value,
    df = denominator_df,
    exact_ems_match = exact_match,
    exact_single_mean_square = exact_single,
    pooled_strata = names(coefficients)[coefficients != 0],
    ems_achieved = setNames(
      as.vector(t(denominator_ems) %*% coefficients), colnames(denominator_ems)
    ),
    ems_target = numerator_ems,
    warnings = if (any(coefficients < -tol)) {
      "Quasi-denominator contains negative mean-square coefficients."
    } else character()
  )
}

.required_fixed_random_interactions <- function(design) {
  if (!length(design$fixed_effects) || !length(design$random_terms)) return(list())
  declared <- vapply(design$random_terms, function(x) {
    .canonical_term(paste(x$variables, collapse = ":"))
  }, character(1))
  setNames(lapply(design$fixed_effects, function(fixed_term) {
    fixed_variables <- strsplit(fixed_term, ":", fixed = TRUE)[[1L]]
    required <- character()
    for (random_term in design$random_terms) {
      if (!length(intersect(fixed_variables, random_term$variables))) {
        required <- c(required, .canonical_term(paste(
          unique(c(fixed_variables, random_term$variables)), collapse = ":"
        )))
      }
    }
    unique(setdiff(required, declared))
  }), design$fixed_effects)
}

.validate_crossed_design <- function(design, assume_additive = FALSE, tol = 1e-10) {
  if (!inherits(design, "permanova_design")) {
    stop("`design` must be produced by `compile_permanova_design()`.", call. = FALSE)
  }
  if (is.null(design$G)) stop("The Design Object must include `distance`.", call. = FALSE)
  if (length(design$random_effects) < 1L) {
    stop("Crossed inference requires at least one random-effect grouping term.", call. = FALSE)
  }
  if (is.null(design$random_terms)) {
    stop("The Design Object lacks random-term provenance; recompile the formula.", call. = FALSE)
  }
  origins <- vapply(design$random_terms, `[[`, character(1), "origin")
  if (any(origins == "nesting")) {
    stop("Phase 3 crossed inference does not support `/`-generated nesting.", call. = FALSE)
  }
  if (length(design$estimability$aliased)) {
    stop("Rank-deficient fixed models are unsupported for crossed inference.", call. = FALSE)
  }

  fixed_variables <- all.vars(design$fixed_formula)
  random_variables <- unique(unlist(lapply(design$random_terms, `[[`, "variables")))
  factor_names <- unique(c(fixed_variables, random_variables))
  missing <- setdiff(factor_names, names(design$data))
  if (length(missing)) stop("Unknown design factors: ", paste(missing, collapse = ", "), call. = FALSE)
  is_factor <- vapply(design$data[factor_names], is.factor, logical(1))
  if (any(!is_factor)) {
    stop("Balanced crossed inference currently requires categorical fixed and random factors.", call. = FALSE)
  }

  cells <- interaction(lapply(design$data[factor_names], droplevels), drop = TRUE)
  counts <- table(cells)
  expected_cells <- prod(vapply(design$data[factor_names], nlevels, integer(1)))
  if (length(counts) != expected_cells || length(unique(as.integer(counts))) != 1L) {
    stop("Crossed inference requires a complete, equally replicated factorial design.", call. = FALSE)
  }
  replicates <- unname(counts[[1L]])
  if (replicates < 2L) {
    stop("At least two observations per crossed factorial cell are required.", call. = FALSE)
  }

  missing_by_term <- .required_fixed_random_interactions(design)
  missing_interactions <- unique(unlist(missing_by_term, use.names = FALSE))
  if (length(missing_interactions) && !isTRUE(assume_additive)) {
    declarations <- paste0("(1 | ", missing_interactions, ")", collapse = " + ")
    stop(
      "Unmodelled fixed-by-random interaction covariance component(s): ",
      paste(missing_interactions, collapse = ", "), ". Declare ", declarations,
      " or set `assume_additive = TRUE` to explicitly assume their variances are zero.",
      call. = FALSE
    )
  }

  full_formula <- as.formula(paste("~", paste(factor_names, collapse = " * ")))
  prepared <- .prepare_fixed_model(full_formula, design$data, tol = tol)
  strata <- setNames(lapply(prepared$cumulative, function(term) {
    list(term = term$term, canonical = .canonical_term(term$term), H = term$H_term,
         df = term$df_term)
  }), vapply(prepared$cumulative, `[[`, character(1), "term"))
  strata$Residual <- list(
    term = "Residual", canonical = "Residual",
    H = diag(design$n) - prepared$full$H,
    df = design$n - prepared$full$rank
  )

  components <- c("residual", design$random_effects)
  covariance_bases <- c(list(residual = diag(design$n)), lapply(design$Z, tcrossprod))
  covariance_bases <- covariance_bases[components]
  for (name in names(strata)) {
    stratum <- strata[[name]]
    stratum$ems <- setNames(vapply(covariance_bases, function(V) {
      trace_product(stratum$H, V) / stratum$df
    }, numeric(1)), components)
    stratum$ss <- trace_hgh(stratum$H, design$G)
    stratum$ms <- stratum$ss / stratum$df
    strata[[name]] <- stratum
  }

  list(
    factor_names = factor_names, fixed_variables = fixed_variables,
    replicates = replicates, strata = strata, components = components,
    missing_by_term = missing_by_term,
    missing_interactions = missing_interactions,
    assume_additive = isTRUE(assume_additive)
  )
}

.validate_crossed_permutation_design <- function(design, structure,
                                                 assume_additive = FALSE,
                                                 tol = 1e-10) {
  if (isTRUE(assume_additive)) {
    stop(
      "Crossed permutation inference does not permit `assume_additive = TRUE`; ",
      "declare the treatment-by-block interaction component.",
      call. = FALSE
    )
  }
  fixed_variables <- all.vars(design$fixed_formula)
  if (length(design$fixed_effects) != 1L || length(fixed_variables) != 1L ||
      .canonical_term(design$fixed_effects[[1L]]) != fixed_variables[[1L]]) {
    stop(
      "Crossed permutation inference supports exactly one categorical fixed treatment factor.",
      call. = FALSE
    )
  }
  treatment_name <- fixed_variables[[1L]]
  random_origins <- vapply(design$random_terms, `[[`, character(1), "origin")
  main_terms <- design$random_terms[random_origins == "main"]
  interaction_terms <- design$random_terms[random_origins == "interaction"]
  if (length(main_terms) != 1L || length(interaction_terms) != 1L ||
      length(design$random_terms) != 2L) {
    stop(
      "Crossed permutation inference requires exactly one random block and its declared treatment-by-block interaction.",
      call. = FALSE
    )
  }
  block_term <- main_terms[[1L]]
  interaction_term <- interaction_terms[[1L]]
  if (length(block_term$variables) != 1L ||
      treatment_name %in% block_term$variables) {
    stop("The random block must be one factor distinct from treatment.",
         call. = FALSE)
  }
  block_name <- block_term$variables[[1L]]
  if (!setequal(interaction_term$variables, c(treatment_name, block_name))) {
    stop(
      "The declared random interaction must contain exactly treatment and block.",
      call. = FALSE
    )
  }
  treatment <- droplevels(design$data[[treatment_name]])
  block <- droplevels(design$data[[block_name]])
  if (!is.factor(treatment) || !is.factor(block) || nlevels(treatment) < 2L ||
      nlevels(block) < 2L) {
    stop("Treatment and block must each be categorical with at least two levels.",
         call. = FALSE)
  }

  cell <- droplevels(design$groups[[interaction_term$label]])
  cell_blocks <- unname(split(seq_len(design$n), cell))
  mapping <- unique(data.frame(
    cell = as.integer(cell), treatment = as.integer(treatment),
    block = as.integer(block)
  ))
  if (any(table(mapping$cell) != 1L) ||
      nrow(mapping) != nlevels(treatment) * nlevels(block)) {
    stop("Every treatment-by-block cell must map to one treatment and one block.",
         call. = FALSE)
  }
  incidence <- table(mapping$block, mapping$treatment)
  if (any(incidence != 1L)) {
    stop("Every treatment must occur exactly once as a cell within every block.",
         call. = FALSE)
  }
  mapping <- mapping[order(mapping$block, mapping$treatment), , drop = FALSE]
  cells_by_block <- unname(split(mapping$cell, mapping$block))
  if (length(unique(lengths(cells_by_block))) != 1L ||
      any(lengths(cells_by_block) != nlevels(treatment))) {
    stop("Each block must contain the same complete set of treatment cells.",
         call. = FALSE)
  }
  if (length(unique(lengths(cell_blocks))) != 1L ||
      lengths(cell_blocks)[[1L]] < 2L) {
    stop("Crossed permutation inference requires equal replication of at least two observations per cell.",
         call. = FALSE)
  }

  canonical_strata <- vapply(structure$strata, `[[`, character(1), "canonical")
  numerator_index <- match(.canonical_term(treatment_name), canonical_strata)
  interaction_index <- match(
    .canonical_term(paste(treatment_name, block_name, sep = ":")),
    canonical_strata
  )
  if (is.na(numerator_index) || is.na(interaction_index)) {
    stop("Treatment or treatment-by-block ANOVA stratum could not be identified.",
         call. = FALSE)
  }
  numerator <- structure$strata[[numerator_index]]
  denominator <- structure$strata[[interaction_index]]
  if (max(abs(numerator$ems - denominator$ems)) > sqrt(tol)) {
    stop("Treatment and treatment-by-block strata do not have matching null EMS.",
         call. = FALSE)
  }

  list(
    treatment_name = treatment_name, block_name = block_name,
    interaction_name = interaction_term$label,
    treatment = treatment, block = block, cell = cell,
    treatments = nlevels(treatment), blocks = nlevels(block),
    replicates = lengths(cell_blocks)[[1L]],
    cell_blocks = cell_blocks, cells_by_block = cells_by_block,
    mapping = mapping, numerator = numerator, denominator = denominator
  )
}

.crossed_block_permutation_plan <- function(permutation_structure,
                                            permutations, max_exact) {
  if (length(max_exact) != 1L || !is.finite(max_exact) || max_exact < 1 ||
      max_exact != as.integer(max_exact)) {
    stop("`max_exact` must be one positive integer.", call. = FALSE)
  }
  attainable <- .safe_factorial_product(
    rep(permutation_structure$treatments, permutation_structure$blocks)
  )
  enumerate <- FALSE
  if (is.character(permutations)) {
    if (length(permutations) != 1L || permutations != "all") {
      stop("`permutations` must be a non-negative integer or \"all\".",
           call. = FALSE)
    }
    if (!is.finite(attainable) || attainable > max_exact) {
      stop("Exact enumeration exceeds `max_exact`; request a finite Monte Carlo sample.",
           call. = FALSE)
    }
    enumerate <- TRUE
  } else {
    if (length(permutations) != 1L || !is.finite(permutations) ||
        permutations < 0 || permutations != as.integer(permutations)) {
      stop("`permutations` must be a non-negative integer or \"all\".",
           call. = FALSE)
    }
    permutations <- as.integer(permutations)
    enumerate <- permutations > 0L && is.finite(attainable) &&
      permutations >= attainable - 1L && attainable <= max_exact
  }

  n_cells <- length(permutation_structure$cell_blocks)
  if (enumerate) {
    local_orders <- lapply(
      permutation_structure$cells_by_block,
      function(cells) .all_orders(length(cells))
    )
    selections <- do.call(expand.grid, c(
      lapply(local_orders, function(x) seq_len(nrow(x))),
      list(KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    ))
    orders <- matrix(
      rep(seq_len(n_cells), each = nrow(selections)),
      nrow = nrow(selections)
    )
    for (block_index in seq_along(local_orders)) {
      cells <- permutation_structure$cells_by_block[[block_index]]
      for (row in seq_len(nrow(selections))) {
        local_order <- local_orders[[block_index]][
          selections[row, block_index], , drop = TRUE
        ]
        orders[row, cells] <- cells[local_order]
      }
    }
  } else if (permutations > 0L) {
    orders <- matrix(rep(seq_len(n_cells), each = permutations),
                     nrow = permutations)
    for (row in seq_len(permutations)) {
      for (cells in permutation_structure$cells_by_block) {
        orders[row, cells] <- sample(cells, length(cells))
      }
    }
  } else {
    orders <- matrix(integer(), nrow = 0L, ncol = n_cells)
  }
  list(
    orders = orders, exact = enumerate, attainable = attainable,
    algorithm = paste0(
      "raw permutation of intact `", permutation_structure$interaction_name,
      "` cells within `", permutation_structure$block_name, "`"
    )
  )
}

.crossed_permutation_statistics <- function(design, permutation_structure,
                                            plan, tol = 1e-10) {
  numerator <- permutation_structure$numerator
  denominator <- permutation_structure$denominator
  observed_denominator <- denominator$ms
  if (.mean_square_is_zero(observed_denominator, design$G, tol)) {
    stop("The treatment-by-block denominator mean square is numerically zero.",
         call. = FALSE)
  }
  observed <- numerator$ms / observed_denominator
  statistics <- .permuted_statistics(
    nrow(plan$orders),
    function(row) {
      index <- .block_permutation_index(
        permutation_structure$cell_blocks, plan$orders[row, ]
      )
      G_permuted <- design$G[index, index, drop = FALSE]
      numerator_ms <- trace_hgh(numerator$H, G_permuted) / numerator$df
      denominator_ms <- trace_hgh(denominator$H, G_permuted) / denominator$df
      numerator_ms / denominator_ms
    },
    label = "Permuting intact treatment cells within blocks"
  )
  inference <- .permutation_inference(statistics, observed, plan$exact, tol)
  p_value <- inference$p_value
  list(observed = observed, statistics = statistics, p_value = p_value)
}

#' Estimate distance-based variance components
#'
#' Uses the residual and each named random-effect ANOVA stratum to solve the
#' balanced trace-derived EMS system. Negative unrestricted estimates are
#' retained; truncation is an optional reported sensitivity analysis.
#'
#' @param design A balanced crossed canonical Design Object.
#' @param truncate Also report estimates truncated at zero.
#' @param assume_additive Explicitly assume omitted fixed-by-random interaction
#'   variances are zero. The default rejects such omissions.
#' @param tol Numerical tolerance.
#' @return A variance-component specification.
#' @export
estimate_variance_components <- function(design, truncate = FALSE,
                                         assume_additive = FALSE, tol = 1e-10) {
  structure <- .validate_crossed_design(
    design, assume_additive = assume_additive, tol = tol
  )
  canonical_names <- vapply(structure$strata, `[[`, character(1), "canonical")
  random_indices <- match(
    vapply(design$random_effects, .canonical_term, character(1)), canonical_names
  )
  if (anyNA(random_indices)) {
    stop("Required random-effect EMS strata could not be identified.", call. = FALSE)
  }
  residual_target <- setNames(c(1, rep(0, length(design$random_effects))), structure$components)
  fixed_indices <- match(
    vapply(design$fixed_effects, .canonical_term, character(1)), canonical_names
  )
  residual_indices <- which(vapply(structure$strata, function(x) {
    max(abs(x$ems - residual_target)) <= tol
  }, logical(1)))
  residual_indices <- setdiff(residual_indices, fixed_indices)
  residual_df <- sum(vapply(structure$strata[residual_indices], `[[`, numeric(1), "df"))
  residual_ms <- sum(vapply(structure$strata[residual_indices], `[[`, numeric(1), "ss")) /
    residual_df
  selected <- structure$strata[random_indices]
  coefficient_matrix <- rbind(
    pooled_residual = residual_target,
    do.call(rbind, lapply(selected, `[[`, "ems"))
  )
  observed_ms <- c(pooled_residual = residual_ms,
                   vapply(selected, `[[`, numeric(1), "ms"))
  estimates <- tryCatch(
    setNames(as.vector(solve(coefficient_matrix, observed_ms)), structure$components),
    error = function(e) stop("Variance-component EMS system is singular.", call. = FALSE)
  )
  notes <- if (any(estimates < -tol)) {
    "Negative unrestricted variance components were retained."
  } else character()
  if (length(structure$missing_interactions)) {
    notes <- c(notes, paste0(
      "Explicit additive assumption: variance is zero for ",
      paste(structure$missing_interactions, collapse = ", "), "."
    ))
  }

  list(
    unrestricted = estimates,
    truncated = if (truncate) pmax(estimates, 0) else NULL,
    coefficient_matrix = coefficient_matrix,
    observed_mean_squares = observed_ms,
    strata = c(
      paste0("pooled: ", paste(vapply(structure$strata[residual_indices], `[[`,
                                      character(1), "term"), collapse = ", ")),
      vapply(selected, `[[`, character(1), "term")
    ),
    warnings = notes,
    assumptions = list(
      additive_opt_in = structure$assume_additive,
      zero_variance_components = structure$missing_interactions
    ),
    method = "balanced trace-derived EMS equations"
  )
}

#' Balanced crossed distance-based mixed model
#'
#' Derives numerical EMS coefficients from covariance and hypothesis
#' projections, chooses exact mean-square denominators where possible, and
#' otherwise constructs an EMS-matched quasi-F denominator.
#'
#' @param design Canonical complete balanced crossed-factor design.
#' @param p_method `"none"` or opt-in `"satterthwaite"` F approximation.
#' @param truncate Also report zero-truncated variance-component estimates.
#' @param assume_additive Explicitly assume omitted fixed-by-random interaction
#'   variances are zero. The default rejects such omissions.
#' @param tol Numerical tolerance.
#' @param permutations A non-negative integer or `"all"` when
#'   `p_method = "permutation"`.
#' @param seed Optional permutation seed.
#' @param keep_permutations Retain permutation statistics.
#' @param max_exact Maximum exact permutation space.
#' @return An object of class `permanova_crossed`.
#' @export
permanova_crossed <- function(design,
                              p_method = c("none", "satterthwaite", "permutation"),
                              truncate = FALSE, assume_additive = FALSE,
                              tol = 1e-10, permutations = 999L, seed = NULL,
                              keep_permutations = FALSE,
                              max_exact = 100000L) {
  p_method <- match.arg(p_method)
  structure <- .validate_crossed_design(
    design, assume_additive = assume_additive, tol = tol
  )
  semimetric <- !is.null(design$spectral) && !isTRUE(design$spectral$euclidean)
  if (p_method == "satterthwaite" && semimetric) {
    warning(
      "Satterthwaite p-values were suppressed because the dissimilarity is non-Euclidean; use a justified permutation test or report a distance-correction sensitivity analysis.",
      call. = FALSE
    )
  }
  permutation_structure <- NULL
  permutation_plan <- NULL
  permutation_result <- NULL
  if (p_method == "permutation") {
    permutation_structure <- .validate_crossed_permutation_design(
      design, structure, assume_additive = assume_additive, tol = tol
    )
    if (!is.null(seed)) set.seed(seed)
    permutation_plan <- .crossed_block_permutation_plan(
      permutation_structure, permutations = permutations,
      max_exact = max_exact
    )
    if (nrow(permutation_plan$orders) > 0L &&
        is.finite(permutation_plan$attainable) &&
        permutation_plan$attainable < 100) {
      warning(
        "Only ", permutation_plan$attainable,
        " within-block treatment-cell permutations are attainable; ",
        "interpret p-values cautiously.", call. = FALSE
      )
    }
    permutation_result <- .crossed_permutation_statistics(
      design, permutation_structure, permutation_plan, tol = tol
    )
  }
  fixed_canonical <- setNames(vapply(design$fixed_effects, .canonical_term, character(1)),
                              design$fixed_effects)
  stratum_canonical <- vapply(structure$strata, `[[`, character(1), "canonical")
  fixed_indices <- match(fixed_canonical, stratum_canonical)
  if (anyNA(fixed_indices)) stop("A fixed-effect hypothesis stratum could not be identified.", call. = FALSE)

  candidate_indices <- setdiff(seq_along(structure$strata), fixed_indices)
  candidate_strata <- structure$strata[candidate_indices]
  candidate_ems <- do.call(rbind, lapply(candidate_strata, `[[`, "ems"))
  rownames(candidate_ems) <- vapply(candidate_strata, `[[`, character(1), "term")
  candidate_ms <- vapply(candidate_strata, `[[`, numeric(1), "ms")
  candidate_df <- vapply(candidate_strata, `[[`, numeric(1), "df")

  tests <- setNames(lapply(seq_along(design$fixed_effects), function(k) {
    term <- design$fixed_effects[[k]]
    numerator <- structure$strata[[fixed_indices[[k]]]]
    denominator <- quasi_denominator(
      numerator$ems, candidate_ems, candidate_ms, candidate_df, tol = tol
    )
    statistic <- numerator$ms / denominator$value
    if (p_method == "permutation" &&
        abs(statistic - permutation_result$observed) > sqrt(tol)) {
      stop("Permutation and EMS engines disagree on the observed pseudo-F.",
           call. = FALSE)
    }
    p_value <- if (p_method == "permutation") {
      permutation_result$p_value
    } else if (p_method == "satterthwaite" && !semimetric && is.finite(statistic) &&
                   statistic >= 0 && is.finite(denominator$df)) {
      stats::pf(statistic, numerator$df, denominator$df, lower.tail = FALSE)
    } else NA_real_
    notes <- denominator$warnings
    missing_for_term <- structure$missing_by_term[[term]]
    assumptions <- character()
    if (length(missing_for_term)) {
      assumptions <- paste0(
        "Explicit additive opt-in assumes zero variance for: ",
        paste(missing_for_term, collapse = ", "), "."
      )
      notes <- c(notes, assumptions)
    }
    if (p_method == "satterthwaite" && semimetric) {
      notes <- c(notes,
        "Satterthwaite p-value suppressed: a non-Euclidean Gower spectrum does not justify the parametric F null.")
    } else if (p_method == "satterthwaite") {
      notes <- c(notes, "P-value uses an approximate parametric Satterthwaite F null distribution.")
    }
    if (p_method == "permutation" && is.finite(permutation_plan$attainable) &&
        permutation_plan$attainable < 100) {
      notes <- c(notes, paste0(
        "Only ", permutation_plan$attainable,
        " within-block treatment-cell permutations are attainable."
      ))
    }
    if (p_method == "permutation" && semimetric) {
      notes <- c(notes,
        "Empirical permutation inference used the retained non-Euclidean Gower spectrum without a parametric F approximation.")
    }
    active <- denominator$coefficients != 0
    expression <- paste0(
      format(denominator$coefficients[active], trim = TRUE), "*MS[",
      names(denominator$coefficients)[active], "]", collapse = " + "
    )
    spec <- list(
      term = term,
      hypothesis = paste0("all estimable coefficients for `", term, "` are zero"),
      numerator = list(ss = numerator$ss, df = numerator$df, ms = numerator$ms,
                       projector = numerator$H),
      denominator = list(
        expression = expression, coefficients = denominator$coefficients,
        ms = denominator$value, df = denominator$df,
        pooled_strata = denominator$pooled_strata,
        assumptions = assumptions
      ),
      symbolic_ems = list(
        components = numerator$ems, denominator_components = denominator$ems_achieved,
        derivation = paste0(
          "EMS coefficients are tr(H Z Z') / df. Matching denominator strata pooled: ",
          paste(denominator$pooled_strata, collapse = ", "), ".",
          if (length(assumptions)) paste0(" ", assumptions) else ""
        ),
        declared_covariance_components = structure$components,
        zero_variance_assumptions = missing_for_term
      ),
      quasi_f_coefficients = denominator$coefficients,
      statistic = statistic, p_value = p_value,
      exchangeable_units = if (p_method == "permutation") {
        permutation_structure$interaction_name
      } else design$random_effects,
      permutation_strata = if (p_method == "permutation") {
        permutation_structure$block_name
      } else NULL,
      permutation_algorithm = if (p_method == "permutation") {
        permutation_plan$algorithm
      } else "none",
      attainable_permutations = if (p_method == "permutation") {
        permutation_plan$attainable
      } else NA_real_,
      permutations_used = if (p_method == "permutation") {
        nrow(permutation_plan$orders)
      } else 0L,
      unique_permutations_used = if (p_method == "permutation" &&
                                      nrow(permutation_plan$orders)) {
        nrow(unique(permutation_plan$orders))
      } else 0L,
      exact_permutation = p_method == "permutation" &&
        permutation_plan$exact && nrow(permutation_plan$orders) > 0L,
      permutation_seed = if (p_method == "permutation") seed else NULL,
      approximation_status = if (p_method == "permutation" &&
                                  !nrow(permutation_plan$orders)) {
        "permutation p-value not requested"
      } else if (p_method == "permutation" && permutation_plan$exact) {
        "exhaustive within-block treatment-cell permutation"
      } else if (p_method == "permutation") {
        "Monte Carlo within-block treatment-cell permutation"
      } else if (p_method == "satterthwaite" && semimetric) {
        "parametric p-value suppressed for non-Euclidean dissimilarity"
      } else if (p_method == "satterthwaite") {
        "Satterthwaite parametric approximation"
      } else if (denominator$exact_ems_match) "exact EMS denominator; p-value not requested" else
        "quasi-F denominator; p-value not requested",
      warnings = notes,
      permutation_statistics = if (p_method == "permutation" &&
                                   keep_permutations) {
        permutation_result$statistics
      } else NULL
    )
    class(spec) <- "permanova_test_spec"
    spec
  }), design$fixed_effects)

  table <- do.call(rbind, lapply(tests, function(spec) data.frame(
    term = spec$term, df = spec$numerator$df, ss = spec$numerator$ss,
    ms = spec$numerator$ms, denominator = spec$denominator$expression,
    denominator_df = spec$denominator$df, pseudo_f = spec$statistic,
    p_value = spec$p_value, stringsAsFactors = FALSE
  )))
  rownames(table) <- NULL

  out <- list(
    call = match.call(), design = design, table = table, tests = tests,
    variance_components = estimate_variance_components(
      design, truncate = truncate, assume_additive = assume_additive, tol = tol
    ),
    strata = structure$strata,
    supported_design = paste0(
      "complete balanced factorial with declared crossed random-intercept interactions",
      if (structure$assume_additive) "; additive omissions explicitly opted in" else ""
    ),
    assumptions = list(
      additive_opt_in = structure$assume_additive,
      zero_variance_components = structure$missing_interactions
    ),
    p_method = p_method,
    permutation = if (p_method == "permutation") list(
      algorithm = permutation_plan$algorithm,
      treatment = permutation_structure$treatment_name,
      block = permutation_structure$block_name,
      cell = permutation_structure$interaction_name,
      treatments = permutation_structure$treatments,
      blocks = permutation_structure$blocks,
      replicates = permutation_structure$replicates,
      exact = permutation_plan$exact,
      attainable = permutation_plan$attainable,
      used = nrow(permutation_plan$orders), seed = seed
    ) else NULL
  )
  class(out) <- "permanova_crossed"
  out
}

#' @export
print.permanova_crossed <- function(x, ...) {
  cat("Distance-based crossed mixed model\n")
  cat("Supported design:", x$supported_design, "\n")
  if (length(x$assumptions$zero_variance_components)) {
    cat("Assumed zero variance:",
        paste(x$assumptions$zero_variance_components, collapse = ", "), "\n")
  }
  print(x$table, row.names = FALSE)
  cat("Unrestricted variance components:\n")
  print(x$variance_components$unrestricted)
  invisible(x)
}
