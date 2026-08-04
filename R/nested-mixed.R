.validate_nested_design <- function(design, tol = 1e-10) {
  if (!inherits(design, "permanova_design")) {
    stop("`design` must be produced by `compile_permanova_design()`.", call. = FALSE)
  }
  if (is.null(design$G)) stop("The Design Object must include `distance`.", call. = FALSE)
  if (length(design$random_terms) < 2L) {
    stop("Nested inference requires a `/`-generated hierarchy with at least two levels.", call. = FALSE)
  }
  origins <- vapply(design$random_terms, `[[`, character(1), "origin")
  if (origins[[1L]] != "main" || any(origins[-1L] != "nesting")) {
    stop("Nested inference supports one pure `/`-generated hierarchy and no crossed random interactions.", call. = FALSE)
  }
  labels <- vapply(design$random_terms, `[[`, character(1), "label")
  parents <- vapply(design$random_terms, function(x) x$parent, character(1))
  if (!all(parents[-1L] == labels[-length(labels)])) {
    stop("Nested random terms must form one contiguous containment chain.", call. = FALSE)
  }
  if (!length(design$fixed_effects)) {
    stop("At least one fixed-effect term is required.", call. = FALSE)
  }
  if (length(design$estimability$aliased) ||
      any(!design$estimability$terms)) {
    stop("The fixed-effect model is rank deficient; aliased nested hypotheses are unsupported.", call. = FALSE)
  }

  top_name <- labels[[1L]]
  top <- droplevels(design$groups[[top_name]])
  top_blocks <- unname(split(seq_len(design$n), top))
  if (length(unique(lengths(top_blocks))) != 1L) {
    stop("Balanced nested inference requires equal observations per top-level unit.", call. = FALSE)
  }
  if (any(!design$balance$random_balanced[labels])) {
    stop("Balanced nested inference requires equal observations in every declared nested unit.", call. = FALSE)
  }

  children_per_parent <- integer(length(labels) - 1L)
  for (level in seq.int(2L, length(labels))) {
    parent <- design$groups[[labels[[level - 1L]]]]
    child <- design$groups[[labels[[level]]]]
    counts <- vapply(split(child, parent), function(x) nlevels(droplevels(x)), integer(1))
    if (length(unique(counts)) != 1L) {
      stop("Balanced nested inference requires equal child-unit counts within every parent.", call. = FALSE)
    }
    children_per_parent[[level - 1L]] <- counts[[1L]]
  }
  deepest_sizes <- table(design$groups[[tail(labels, 1L)]])
  if (min(deepest_sizes) < 2L) {
    stop("At least two observations per deepest nested unit are required.", call. = FALSE)
  }

  components <- c("residual", labels)
  covariance_bases <- c(list(residual = diag(design$n)), lapply(design$Z[labels], tcrossprod))
  covariance_bases <- covariance_bases[components]
  ems_for <- function(H, df) setNames(vapply(covariance_bases, function(V) {
    trace_product(H, V) / df
  }, numeric(1)), components)

  intercept_projection <- projection_matrix(
    matrix(1, nrow = design$n, ncol = 1L), tol = tol
  )
  boundary_projections <- c(
    list(intercept_projection), design$projections$random[labels],
    list(list(H = diag(design$n), rank = design$n))
  )
  level_names <- c(labels, "observation")
  parent_names <- c("grand mean", labels)
  raw_strata <- vector("list", length(level_names))
  for (level in seq_along(level_names)) {
    parent_projection <- boundary_projections[[level]]
    unit_projection <- boundary_projections[[level + 1L]]
    nesting_error <- max(abs(
      unit_projection$H %*% parent_projection$H - parent_projection$H
    ))
    if (nesting_error > max(100 * tol, sqrt(.Machine$double.eps))) {
      stop("Declared nested grouping projections do not form a containment chain.",
           call. = FALSE)
    }
    H_raw <- unit_projection$H - parent_projection$H
    raw_df <- unit_projection$rank - parent_projection$rank
    projected_X <- H_raw %*% design$X
    projected_norm <- sqrt(colSums(projected_X ^ 2))
    original_norm <- pmax(sqrt(colSums(design$X ^ 2)), 1)
    projected_X[, projected_norm / original_norm <=
                  max(100 * tol, sqrt(.Machine$double.eps))] <- 0
    fixed_within <- projection_matrix(projected_X, tol = tol)
    denominator_df <- raw_df - fixed_within$rank
    if (denominator_df <= 0L) {
      stop(
        "No residual degrees of freedom remain in hierarchy stratum `",
        level_names[[level]], " within ", parent_names[[level]], "`.",
        call. = FALSE
      )
    }
    H_denominator <- H_raw - fixed_within$H
    projector_error <- max(abs(
      H_denominator %*% H_denominator - H_denominator
    ))
    if (projector_error > max(100 * tol, sqrt(.Machine$double.eps))) {
      stop(
        "Fixed effects do not decompose cleanly across the balanced nested strata.",
        call. = FALSE
      )
    }
    raw_strata[[level]] <- list(
      level = level, level_name = level_names[[level]],
      parent_name = parent_names[[level]], H_raw = H_raw, raw_df = raw_df,
      H_fixed = fixed_within$H, fixed_df = fixed_within$rank,
      H = H_denominator, df = denominator_df
    )
  }

  assignment_tolerance <- max(100 * tol, sqrt(.Machine$double.eps))
  term_levels <- setNames(integer(length(design$fixed_effects)), design$fixed_effects)
  term_assignment_error <- setNames(numeric(length(design$fixed_effects)),
                                    design$fixed_effects)
  for (term in design$fixed_effects) {
    H_term <- design$projections$terms[[term]]$H_term
    errors <- vapply(raw_strata, function(stratum) {
      max(abs(H_term - stratum$H_raw %*% H_term %*% stratum$H_raw))
    }, numeric(1))
    matches <- which(errors <= assignment_tolerance)
    if (length(matches) != 1L) {
      fixed_variables <- all.vars(design$fixed_formula)
      continuous <- fixed_variables[vapply(fixed_variables, function(variable) {
        value <- design$data[[variable]]
        is.numeric(value) && !is.factor(value)
      }, logical(1))]
      advice <- if (length(continuous)) {
        paste0(
          "The fixed model includes continuous covariate(s) (",
          paste(continuous, collapse = ", "),
          ") varying within the nested units. Covariate adjustment in a nested ",
          "or mixed design is not supported by `permanova_nested()`; the ",
          "balanced-EMS engine requires every fixed term to occupy one hierarchy ",
          "stratum. Remove the covariate, or reduce it to a single-level ",
          "predictor."
        )
      } else {
        paste0(
          "Use a level-specific, within-parent-balanced parameterization. For ",
          "balanced cross-level factor interactions, assign explicit orthogonal ",
          "contrasts such as `contr.sum` to the crossed factors before compiling ",
          "the design."
        )
      }
      stop(
        "Fixed term `", term,
        "` spans multiple hierarchy levels or is not orthogonal to its parent ",
        "stratum. ", advice, call. = FALSE
      )
    }
    term_levels[[term]] <- matches[[1L]]
    term_assignment_error[[term]] <- errors[[matches[[1L]]]]
  }

  strata <- list()
  for (level in seq_along(raw_strata)) {
    stratum <- raw_strata[[level]]
    name <- if (level == 1L) {
      paste0(level_names[[level]], " within fixed model")
    } else if (level <= length(labels)) {
      paste0(level_names[[level]], " within ", parent_names[[level]],
             " and fixed model")
    } else {
      paste0("Within ", tail(labels, 1L), " after fixed model")
    }
    stratum$term <- name
    stratum$assigned_terms <- names(term_levels)[term_levels == level]
    stratum$ss <- trace_hgh(stratum$H, design$G)
    stratum$ms <- stratum$ss / stratum$df
    stratum$ems <- ems_for(stratum$H, stratum$df)
    strata[[level]] <- stratum
  }
  names(strata) <- vapply(strata, `[[`, character(1), "term")

  unit_factors <- c(
    design$groups[labels],
    list(observation = factor(seq_len(design$n)))
  )

  list(
    labels = labels, top_name = top_name, top = top, top_blocks = top_blocks,
    top_units = length(top_blocks), observations_per_top = lengths(top_blocks)[[1L]],
    children_per_parent = setNames(children_per_parent, labels[-1L]),
    deepest_replicates = unname(deepest_sizes[[1L]]),
    level_names = level_names, parent_names = parent_names,
    unit_factors = unit_factors, raw_strata = raw_strata,
    term_levels = term_levels,
    term_assignment_error = term_assignment_error,
    H_denominator = strata[[1L]]$H, denominator_df = strata[[1L]]$df,
    denominator_name = names(strata)[[1L]],
    components = components, covariance_bases = covariance_bases,
    ems_for = ems_for, strata = strata
  )
}

.safe_factorial_product <- function(counts) {
  log_value <- sum(lgamma(counts + 1))
  if (!is.finite(log_value) || log_value > log(.Machine$double.xmax)) {
    Inf
  } else {
    prod(factorial(counts))
  }
}

.nested_exchangeability <- function(structure, level) {
  unit <- droplevels(structure$unit_factors[[level]])
  unit_blocks <- unname(split(seq_along(unit), unit))
  if (level == 1L) {
    units_by_parent <- list(seq_along(unit_blocks))
  } else {
    parent <- droplevels(structure$unit_factors[[level - 1L]])
    mapping <- unique(data.frame(
      unit = as.integer(unit), parent = as.integer(parent)
    ))
    if (any(table(mapping$unit) != 1L)) {
      stop("A nested exchangeability unit belongs to more than one parent.",
           call. = FALSE)
    }
    units_by_parent <- unname(split(mapping$unit, mapping$parent))
  }
  counts <- lengths(units_by_parent)
  if (length(unique(counts)) != 1L) {
    stop("Balanced nested permutations require equal exchangeable units per parent.",
         call. = FALSE)
  }
  list(
    level = level, level_name = structure$level_names[[level]],
    parent_name = structure$parent_names[[level]],
    unit_blocks = unit_blocks, units_by_parent = units_by_parent,
    units_per_parent = counts[[1L]], parents = length(units_by_parent),
    attainable = .safe_factorial_product(counts)
  )
}

.nested_permutation_plan <- function(structure, level, permutations,
                                     max_exact) {
  exchangeability <- .nested_exchangeability(structure, level)
  attainable <- exchangeability$attainable
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

  n_units <- length(exchangeability$unit_blocks)
  if (enumerate) {
    local_orders <- lapply(
      exchangeability$units_by_parent,
      function(units) .all_orders(length(units))
    )
    selections <- do.call(expand.grid, c(
      lapply(local_orders, function(x) seq_len(nrow(x))),
      list(KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    ))
    orders <- matrix(
      rep(seq_len(n_units), each = nrow(selections)),
      nrow = nrow(selections)
    )
    for (parent_index in seq_along(local_orders)) {
      units <- exchangeability$units_by_parent[[parent_index]]
      for (row in seq_len(nrow(selections))) {
        local_order <- local_orders[[parent_index]][
          selections[row, parent_index], , drop = TRUE
        ]
        orders[row, units] <- units[local_order]
      }
    }
  } else if (permutations > 0L) {
    orders <- matrix(rep(seq_len(n_units), each = permutations),
                     nrow = permutations)
    for (row in seq_len(permutations)) {
      for (units in exchangeability$units_by_parent) {
        orders[row, units] <- sample(units, length(units))
      }
    }
  } else {
    orders <- matrix(integer(), nrow = 0L, ncol = n_units)
  }

  c(exchangeability, list(
    orders = orders, exact = enumerate,
    algorithm = paste0(
      "Freedman-Lane residual permutation of `",
      exchangeability$level_name, "` units within `",
      exchangeability$parent_name, "`, moving descendants intact"
    )
  ))
}

#' Derive EMS for a balanced nested random-intercept hierarchy
#'
#' @param design Canonical Design Object with a `/`-generated hierarchy.
#' @param tol Numerical tolerance.
#' @return A named list of exact nested EMS specifications.
#' @export
derive_nested_ems <- function(design, tol = 1e-10) {
  structure <- .validate_nested_design(design, tol = tol)
  setNames(lapply(design$fixed_effects, function(term) {
    projection <- design$projections$terms[[term]]
    level <- structure$term_levels[[term]]
    denominator_stratum <- structure$strata[[level]]
    denominator_ems <- denominator_stratum$ems
    numerator_ems <- structure$ems_for(projection$H_term, projection$df)
    mismatch <- max(abs(numerator_ems - denominator_ems))
    if (mismatch > sqrt(tol)) {
      stop("No exact nested denominator matches fixed term `", term, "`.", call. = FALSE)
    }
    list(
      term = term,
      hypothesis = paste0("all estimable coefficients for `", term, "` are zero"),
      numerator_df = projection$df,
      numerator_under_null = numerator_ems,
      numerator_under_alternative = c(numerator_ems, fixed_effect = 1),
      denominator = denominator_ems,
      denominator_term = denominator_stratum$term,
      denominator_df = denominator_stratum$df,
      assignment_level = denominator_stratum$level_name,
      assignment_level_index = level,
      parent_level = denominator_stratum$parent_name,
      assignment_error = structure$term_assignment_error[[term]],
      exact_denominator = TRUE,
      derivation = paste0(
        "Term `", term, "` lies wholly in the `",
        denominator_stratum$level_name, " within ",
        denominator_stratum$parent_name, "` hierarchy stratum. Its numerator and `",
        denominator_stratum$term,
        "` denominator have identical trace-derived residual and nested ",
        "random-intercept EMS coefficients under H0."
      )
    )
  }), design$fixed_effects)
}

.nested_test_specification <- function(term, design, structure, ems, plan,
                                       tol) {
  projection <- design$projections$terms[[term]]
  level <- structure$term_levels[[term]]
  denominator_stratum <- structure$strata[[level]]
  numerator_ss <- trace_hgh(projection$H_term, design$G)
  numerator_ms <- numerator_ss / projection$df
  denominator_ms <- denominator_stratum$ms
  if (abs(denominator_ms) <= tol) {
    stop("The nested denominator mean square is numerically zero for term `", term, "`.", call. = FALSE)
  }
  observed <- numerator_ms / denominator_ms
  permuted <- .freedman_lane_statistics(
    design$G, projection$H_reduced,
    projection$H_term, projection$df,
    denominator_stratum$H, denominator_stratum$df,
    plan$orders,
    index = function(order) .block_permutation_index(plan$unit_blocks, order),
    label = paste0("Permuting term `", term, "`"), tol = tol
  )
  p_value <- NA_real_
  if (length(permuted)) {
    extreme <- sum(permuted >= observed - tol, na.rm = TRUE)
    p_value <- if (plan$exact) extreme / length(permuted) else
      (extreme + 1) / (length(permuted) + 1)
  }
  notes <- character()
  if (is.finite(plan$attainable) && plan$attainable < 100) {
    notes <- paste0(
      "Only ", plan$attainable, " level-specific permutations are attainable."
    )
  }
  if (numerator_ss < -tol || denominator_stratum$ss < -tol) {
    notes <- c(notes, "A negative distance-based sum of squares was retained.")
  }

  out <- list(
    term = term, hypothesis = ems$hypothesis,
    numerator = list(
      ss = numerator_ss, df = projection$df, ms = numerator_ms,
      projector = projection$H_term
    ),
    denominator = list(
      expression = paste0("MS[", denominator_stratum$term, "]"),
      coefficients = 1, ss = denominator_stratum$ss,
      df = denominator_stratum$df, ms = denominator_ms,
      projector = denominator_stratum$H
    ),
    symbolic_ems = ems, quasi_f_coefficients = NULL,
    statistic = observed, p_value = p_value,
    assignment_level = denominator_stratum$level_name,
    assignment_level_index = level,
    parent_level = denominator_stratum$parent_name,
    exchangeable_units = plan$level_name,
    permutation_algorithm = plan$algorithm,
    attainable_permutations = plan$attainable,
    permutations_used = nrow(plan$orders),
    unique_permutations_used = if (nrow(plan$orders)) {
      nrow(unique(plan$orders))
    } else 0L,
    exact_permutation = plan$exact && nrow(plan$orders) > 0L,
    approximation_status = if (!nrow(plan$orders)) {
      "not tested"
    } else if (plan$exact) {
      "exhaustive level-specific permutation"
    } else "Monte Carlo p-value",
    warnings = notes, permutation_statistics = permuted
  )
  class(out) <- "permanova_test_spec"
  out
}

#' Estimate balanced nested distance-based variance components
#'
#' @param design Canonical nested Design Object.
#' @param truncate Also return a separate zero-truncated sensitivity estimate.
#' @param tol Numerical tolerance.
#' @return A nested variance-component specification.
#' @export
estimate_nested_variance_components <- function(design, truncate = FALSE,
                                                tol = 1e-10) {
  structure <- .validate_nested_design(design, tol = tol)
  selected_names <- rev(names(structure$strata))
  selected <- structure$strata[selected_names]
  coefficient_matrix <- do.call(rbind, lapply(selected, `[[`, "ems"))
  observed_ms <- vapply(selected, `[[`, numeric(1), "ms")
  if (nrow(coefficient_matrix) != ncol(coefficient_matrix) ||
      qr(coefficient_matrix, tol = tol)$rank < ncol(coefficient_matrix)) {
    stop("Nested variance-component EMS system is singular.", call. = FALSE)
  }
  estimates <- setNames(
    as.vector(solve(coefficient_matrix, observed_ms)), structure$components
  )
  notes <- if (any(estimates < -tol)) {
    "Negative unrestricted nested variance components were retained."
  } else character()
  list(
    unrestricted = estimates,
    truncated = if (truncate) pmax(estimates, 0) else NULL,
    coefficient_matrix = coefficient_matrix,
    observed_mean_squares = observed_ms,
    strata = selected_names,
    warnings = notes,
    method = "balanced trace-derived nested EMS equations"
  )
}

#' Balanced nested distance-based mixed model
#'
#' Assigns each fixed term to one orthogonal hierarchy stratum, uses that
#' stratum's exact residual mean square, and permutes units within their parent
#' while moving all descendants intact.
#'
#' @param design Canonical Design Object compiled with `distance`.
#' @param permutations A non-negative integer or `"all"`.
#' @param seed Optional random seed.
#' @param keep_permutations Retain permutation statistics.
#' @param truncate Also report zero-truncated nested variance components.
#' @param max_exact Maximum exact permutation space.
#' @param tol Numerical tolerance.
#' @return An object of class `permanova_nested`.
#' @export
permanova_nested <- function(design, permutations = 999L, seed = NULL,
                             keep_permutations = FALSE, truncate = FALSE,
                             max_exact = 100000L, tol = 1e-10) {
  structure <- .validate_nested_design(design, tol = tol)
  ems <- derive_nested_ems(design, tol = tol)
  if (!is.null(seed)) set.seed(seed)
  tested_levels <- unique(unname(structure$term_levels))
  level_plans <- setNames(lapply(tested_levels, function(level) {
    .nested_permutation_plan(
      structure, level, permutations, max_exact
    )
  }), as.character(tested_levels))
  plans <- setNames(lapply(design$fixed_effects, function(term) {
    level_plans[[as.character(structure$term_levels[[term]])]]
  }), design$fixed_effects)
  for (level_key in names(level_plans)) {
    plan <- level_plans[[level_key]]
    if (nrow(plan$orders) > 0L && is.finite(plan$attainable) &&
        plan$attainable < 100) {
      affected_terms <- names(structure$term_levels)[
        structure$term_levels == as.integer(level_key)
      ]
      warning(
        "Only ", plan$attainable,
        " permutations are attainable for hierarchy level `",
        plan$level_name, "` (term(s): ",
        paste(affected_terms, collapse = ", "),
        "); interpret p-values cautiously.", call. = FALSE
      )
    }
  }
  tests <- setNames(lapply(design$fixed_effects, function(term) {
    .nested_test_specification(
      term, design, structure, ems[[term]], plans[[term]], tol
    )
  }), design$fixed_effects)
  if (!keep_permutations) {
    tests <- lapply(tests, function(spec) {
      spec$permutation_statistics <- NULL
      class(spec) <- "permanova_test_spec"
      spec
    })
  }

  table <- do.call(rbind, lapply(tests, function(spec) data.frame(
    term = spec$term, df = spec$numerator$df, ss = spec$numerator$ss,
    ms = spec$numerator$ms, denominator = spec$denominator$expression,
    denominator_df = spec$denominator$df, pseudo_f = spec$statistic,
    p_perm = spec$p_value, assignment_level = spec$assignment_level,
    stringsAsFactors = FALSE
  )))
  for (stratum in structure$strata) {
    table <- rbind(table, data.frame(
      term = stratum$term, df = stratum$df, ss = stratum$ss, ms = stratum$ms,
      denominator = NA_character_, denominator_df = NA_real_,
      pseudo_f = NA_real_, p_perm = NA_real_,
      assignment_level = stratum$level_name
    ))
  }
  rownames(table) <- NULL
  permutation_by_term <- lapply(plans, function(plan) list(
    algorithm = plan$algorithm, assignment_level = plan$level_name,
    parent_level = plan$parent_name, exact = plan$exact,
    attainable = plan$attainable, used = nrow(plan$orders)
  ))
  one_term <- length(permutation_by_term) == 1L
  out <- list(
    call = match.call(), design = design, table = table, tests = tests, ems = ems,
    variance_components = estimate_nested_variance_components(
      design, truncate = truncate, tol = tol
    ),
    hierarchy = structure$labels,
    assignment_levels = setNames(
      structure$level_names[structure$term_levels], names(structure$term_levels)
    ),
    supported_design = paste(
      "balanced nested random intercepts; every fixed term confined to one",
      "orthogonal hierarchy stratum"
    ),
    permutation = list(
      algorithm = "term-specific Freedman-Lane permutation at the assigned hierarchy level",
      exact = if (one_term) permutation_by_term[[1L]]$exact else
        vapply(permutation_by_term, `[[`, logical(1), "exact"),
      attainable = if (one_term) permutation_by_term[[1L]]$attainable else
        vapply(permutation_by_term, `[[`, numeric(1), "attainable"),
      used = if (one_term) permutation_by_term[[1L]]$used else
        vapply(permutation_by_term, `[[`, integer(1), "used"),
      by_term = permutation_by_term, seed = seed
    )
  )
  class(out) <- "permanova_nested"
  out
}

#' @export
print.permanova_nested <- function(x, ...) {
  cat("Distance-based balanced nested mixed model\n")
  cat("Hierarchy:", paste(x$hierarchy, collapse = " > "), "\n")
  cat("Permutation:", x$permutation$algorithm, "\n")
  print(x$table, row.names = FALSE)
  cat("Unrestricted variance components:\n")
  print(x$variance_components$unrestricted)
  invisible(x)
}
