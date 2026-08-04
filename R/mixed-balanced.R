.validate_balanced_random_intercept <- function(design, tol = 1e-10) {
  if (!inherits(design, "permanova_design")) {
    stop("`design` must be produced by `compile_permanova_design()`.", call. = FALSE)
  }
  if (is.null(design$G)) {
    stop("The Design Object must be compiled with `distance` for inference.", call. = FALSE)
  }
  if (length(design$random_effects) != 1L) {
    stop(
      "Phase 2 supports exactly one random-intercept grouping factor; nested or crossed random effects are not yet supported.",
      call. = FALSE
    )
  }
  if (length(design$fixed_effects) == 0L) {
    stop("At least one fixed-effect term is required for testing.", call. = FALSE)
  }
  if (length(design$estimability$aliased)) {
    stop("The fixed-effect model is rank deficient; aliased hypotheses are unsupported.", call. = FALSE)
  }

  cluster_name <- design$random_effects[[1L]]
  cluster <- droplevels(design$groups[[cluster_name]])
  blocks <- split(seq_len(design$n), cluster)
  sizes <- lengths(blocks)
  if (length(unique(sizes)) != 1L) {
    stop("Balanced random-intercept inference requires equal observations per grouping unit.", call. = FALSE)
  }
  if (sizes[[1L]] < 2L) {
    stop("At least two observations per grouping unit are required to identify within-unit variation.", call. = FALSE)
  }

  varies_within <- vapply(seq_len(ncol(design$X)), function(column) {
    any(vapply(blocks, function(index) {
      diff(range(design$X[index, column])) > tol
    }, logical(1)))
  }, logical(1))
  if (any(varies_within)) {
    stop(
      "All fixed predictors must be constant within the random-intercept grouping units for the Phase 2 whole-plot test.",
      call. = FALSE
    )
  }

  H_cluster <- design$projections$random[[cluster_name]]$H
  H_fixed <- design$projections$fixed$H
  rank_cluster <- design$projections$random[[cluster_name]]$rank
  rank_fixed <- design$projections$fixed$rank
  denominator_df <- rank_cluster - rank_fixed
  if (denominator_df <= 0L) {
    stop("No between-grouping-unit residual degrees of freedom remain for an exact denominator.", call. = FALSE)
  }

  list(
    cluster_name = cluster_name,
    cluster = cluster,
    blocks = unname(blocks),
    groups = length(blocks),
    replicates = unname(sizes[[1L]]),
    H_cluster = H_cluster,
    H_denominator = H_cluster - H_fixed,
    denominator_df = denominator_df,
    H_within = diag(design$n) - H_cluster,
    within_df = design$n - rank_cluster
  )
}

.mixed_term_projections <- function(design, hypothesis = "I", tol = 1e-10) {
  hypothesis <- match.arg(hypothesis, c("I", "II", "III"))
  if (identical(hypothesis, "III")) {
    .validate_type3_contrasts(design, tol = tol)
  }
  X <- design$X
  assignment <- attr(X, "assign")
  labels <- design$fixed_effects
  intercept <- which(assignment == 0L)
  setNames(lapply(seq_along(labels), function(index) {
    adjust <- .ss_type_adjust(index, labels, hypothesis)
    reduced_columns <- c(intercept, which(assignment %in% adjust))
    term_columns <- which(assignment == index)
    reduced <- projection_matrix(X[, reduced_columns, drop = FALSE], tol = tol)
    full <- projection_matrix(
      X[, c(reduced_columns, term_columns), drop = FALSE], tol = tol
    )
    df <- full$rank - reduced$rank
    if (df <= 0L) {
      stop("Fixed term `", labels[[index]], "` is not estimable for Type ",
           hypothesis, " mixed-model inference.", call. = FALSE)
    }
    list(
      term = labels[[index]], H_reduced = reduced$H,
      H_term = full$H - reduced$H, df = df, estimable = TRUE
    )
  }), labels)
}

.safe_factorial <- function(n) {
  if (n > 170L) Inf else factorial(n)
}

.all_orders <- function(n) {
  recurse <- function(prefix, remaining) {
    if (!length(remaining)) return(matrix(prefix, nrow = 1L))
    do.call(rbind, lapply(seq_along(remaining), function(k) {
      recurse(c(prefix, remaining[[k]]), remaining[-k])
    }))
  }
  result <- recurse(integer(), seq_len(n))
  matrix(result, ncol = n)
}

.block_permutation_index <- function(blocks, order) {
  index <- seq_len(sum(lengths(blocks)))
  for (target in seq_along(blocks)) index[blocks[[target]]] <- blocks[[order[[target]]]]
  index
}

.permutation_orders <- function(groups, permutations, max_exact) {
  attainable <- .safe_factorial(groups)
  if (is.character(permutations)) {
    if (length(permutations) != 1L || permutations != "all") {
      stop("`permutations` must be a non-negative integer or \"all\".", call. = FALSE)
    }
    if (!is.finite(attainable) || attainable > max_exact) {
      stop("Exact enumeration exceeds `max_exact`; request a finite Monte Carlo sample.", call. = FALSE)
    }
    return(list(orders = .all_orders(groups), exact = TRUE, attainable = attainable))
  }
  if (length(permutations) != 1L || !is.finite(permutations) ||
      permutations < 0 || permutations != as.integer(permutations)) {
    stop("`permutations` must be a non-negative integer or \"all\".", call. = FALSE)
  }

  permutations <- as.integer(permutations)
  if (permutations > 0L && is.finite(attainable) &&
      permutations >= attainable - 1L && attainable <= max_exact) {
    return(list(orders = .all_orders(groups), exact = TRUE, attainable = attainable))
  }
  orders <- if (permutations) {
    t(replicate(permutations, sample.int(groups), simplify = TRUE))
  } else {
    matrix(integer(), nrow = 0L, ncol = groups)
  }
  list(orders = orders, exact = FALSE, attainable = attainable)
}

#' Derive EMS for a balanced random-intercept design
#'
#' The supported Phase 2 envelope contains exactly one random intercept,
#' equal replication per grouping unit, and fixed predictors that are constant
#' within grouping units.
#'
#' @param design A canonical `permanova_design` containing a distance matrix.
#' @param hypothesis Fixed-effect sum-of-squares partition: `"I"`, `"II"`, or
#'   `"III"`.
#' @param tol Numerical tolerance.
#' @return A named list of symbolic EMS specifications.
#' @export
derive_ems <- function(design, hypothesis = "I", tol = 1e-10) {
  structure <- .validate_balanced_random_intercept(design, tol = tol)
  projections <- .mixed_term_projections(design, hypothesis, tol = tol)
  component <- c(residual = 1, random_intercept = structure$replicates)
  setNames(lapply(design$fixed_effects, function(term) {
    projection <- projections[[term]]
    list(
      term = term,
      hypothesis = paste0("all estimable coefficients for `", term, "` are zero"),
      numerator_df = projection$df,
      numerator_under_null = component,
      numerator_under_alternative = c(component, fixed_effect = 1),
      denominator = component,
      denominator_term = paste0(structure$cluster_name, " within fixed model"),
      denominator_df = structure$denominator_df,
      exact_denominator = TRUE,
      derivation = paste0(
        "Both numerator and denominator contain sigma[e]^2 + ",
        structure$replicates, " sigma[", structure$cluster_name,
        "]^2 under H0; the fixed-effect component occurs only in the numerator under H1."
      )
    )
  }), design$fixed_effects)
}

.make_test_specification <- function(term, projection, design, structure, ems, G,
                                     orders, exact, attainable, tol) {
  numerator_ss <- trace_hgh(projection$H_term, G)
  denominator_ss <- trace_hgh(structure$H_denominator, G)
  numerator_ms <- numerator_ss / projection$df
  denominator_ms <- denominator_ss / structure$denominator_df
  if (abs(denominator_ms) <= tol) {
    stop("The exact denominator mean square is numerically zero for term `", term, "`.", call. = FALSE)
  }
  observed <- numerator_ms / denominator_ms

  permuted <- .freedman_lane_statistics(
    G, projection$H_reduced,
    projection$H_term, projection$df,
    structure$H_denominator, structure$denominator_df,
    orders,
    index = function(order) .block_permutation_index(structure$blocks, order),
    label = paste0("Permuting whole plots for term `", term, "`"), tol = tol
  )

  p_value <- NA_real_
  if (length(permuted)) {
    extreme <- sum(permuted >= observed - tol, na.rm = TRUE)
    p_value <- if (exact) extreme / length(permuted) else (extreme + 1) / (length(permuted) + 1)
  }
  notes <- character()
  if (is.finite(attainable) && attainable < 100) {
    notes <- c(notes, paste0("Only ", attainable, " whole-plot permutations are attainable."))
  }
  if (numerator_ss < -tol || denominator_ss < -tol) {
    notes <- c(notes, "A negative distance-based sum of squares was retained.")
  }

  out <- list(
    term = term,
    hypothesis = ems$hypothesis,
    numerator = list(ss = numerator_ss, df = projection$df, ms = numerator_ms,
                     projector = projection$H_term),
    denominator = list(
      expression = paste0("MS[", structure$cluster_name, " within fixed model]"),
      coefficients = 1, ss = denominator_ss, df = structure$denominator_df,
      ms = denominator_ms, projector = structure$H_denominator
    ),
    symbolic_ems = ems,
    quasi_f_coefficients = NULL,
    statistic = observed,
    p_value = p_value,
    exchangeable_units = structure$cluster_name,
    permutation_algorithm = "Freedman-Lane residual permutation of whole plots",
    attainable_permutations = attainable,
    permutations_used = nrow(orders),
    unique_permutations_used = if (nrow(orders)) nrow(unique(orders)) else 0L,
    exact_permutation = exact && nrow(orders) > 0L,
    approximation_status = if (!nrow(orders)) "not tested" else if (exact) "exact" else "Monte Carlo p-value",
    warnings = notes,
    permutation_statistics = permuted
  )
  class(out) <- "permanova_test_spec"
  out
}

#' Balanced distance-based random-intercept inference
#'
#' Builds an explicit test specification for every fixed term, selects the
#' exact between-whole-plot residual denominator from its EMS, and applies
#' hypothesis-specific whole-plot residual permutations.
#'
#' @param design A canonical `permanova_design` compiled with `distance`.
#' @param permutations A non-negative integer or `"all"`. Enumeration is used
#'   automatically when the requested sample covers the attainable space.
#' @param seed Optional random seed.
#' @param hypothesis Fixed-effect sum-of-squares partition: `"I"` (sequential),
#'   `"II"` (marginal), or `"III"` (partial).
#' @param keep_permutations Retain permutation statistics in each test spec.
#' @param max_exact Maximum permutation space allowed for enumeration.
#' @param tol Numerical tolerance.
#' @return An object of class `permanova_mixed`.
#' @export
permanova_mixed <- function(design, permutations = 999L, seed = NULL,
                            hypothesis = "I",
                            keep_permutations = FALSE, max_exact = 100000L,
                            tol = 1e-10) {
  hypothesis <- match.arg(hypothesis, c("I", "II", "III"))
  structure <- .validate_balanced_random_intercept(design, tol = tol)
  projections <- .mixed_term_projections(design, hypothesis, tol = tol)
  ems <- derive_ems(design, hypothesis = hypothesis, tol = tol)
  if (!is.null(seed)) set.seed(seed)
  permutation_plan <- .permutation_orders(structure$groups, permutations, max_exact)
  if (is.finite(permutation_plan$attainable) && permutation_plan$attainable < 100) {
    warning(
      "Only ", permutation_plan$attainable,
      " whole-plot permutations are attainable; interpret p-values cautiously.",
      call. = FALSE
    )
  }

  specifications <- setNames(lapply(design$fixed_effects, function(term) {
    .make_test_specification(
      term, projections[[term]], design, structure, ems[[term]], design$G,
      permutation_plan$orders, permutation_plan$exact,
      permutation_plan$attainable, tol
    )
  }), design$fixed_effects)
  if (!keep_permutations) {
    specifications <- lapply(specifications, function(spec) {
      spec$permutation_statistics <- NULL
      class(spec) <- "permanova_test_spec"
      spec
    })
  }

  table <- do.call(rbind, lapply(specifications, function(spec) {
    data.frame(
      term = spec$term, df = spec$numerator$df, ss = spec$numerator$ss,
      ms = spec$numerator$ms, denominator = spec$denominator$expression,
      denominator_df = spec$denominator$df, pseudo_f = spec$statistic,
      p_perm = spec$p_value, stringsAsFactors = FALSE
    )
  }))
  table <- rbind(table, data.frame(
    term = paste0(structure$cluster_name, " within fixed model"),
    df = structure$denominator_df,
    ss = trace_hgh(structure$H_denominator, design$G),
    ms = trace_hgh(structure$H_denominator, design$G) / structure$denominator_df,
    denominator = NA_character_, denominator_df = NA_integer_,
    pseudo_f = NA_real_, p_perm = NA_real_
  ), data.frame(
    term = paste0("Within ", structure$cluster_name), df = structure$within_df,
    ss = trace_hgh(structure$H_within, design$G),
    ms = trace_hgh(structure$H_within, design$G) / structure$within_df,
    denominator = NA_character_, denominator_df = NA_integer_,
    pseudo_f = NA_real_, p_perm = NA_real_
  ))
  rownames(table) <- NULL

  out <- list(
    call = match.call(), design = design, table = table,
    tests = specifications, ems = ems,
    supported_design = "balanced single random intercept; fixed predictors constant within whole plots",
    permutation = list(
      algorithm = "term-specific Freedman-Lane whole-plot residual permutation",
      exact = permutation_plan$exact,
      attainable = permutation_plan$attainable,
      used = nrow(permutation_plan$orders), seed = seed
    )
  )
  class(out) <- "permanova_mixed"
  out
}

#' @export
print.permanova_mixed <- function(x, ...) {
  cat("Distance-based balanced random-intercept model\n")
  cat("Supported design:", x$supported_design, "\n")
  cat("Permutation:", x$permutation$algorithm, "\n")
  print(x$table, row.names = FALSE)
  invisible(x)
}
