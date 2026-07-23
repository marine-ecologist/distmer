.sign_flip_plan <- function(units, permutations, max_exact) {
  if (length(units) != 1L || !is.finite(units) || units < 1L ||
      units != as.integer(units)) {
    stop("`units` must be one positive integer.", call. = FALSE)
  }
  units <- as.integer(units)
  attainable <- if (units > 1023L) Inf else 2 ^ units
  exact_requested <- is.character(permutations)
  if (exact_requested &&
      (length(permutations) != 1L || permutations != "all")) {
    stop("`permutations` must be a non-negative integer or \"all\".",
         call. = FALSE)
  }
  if (exact_requested && (!is.finite(attainable) || attainable > max_exact)) {
    stop("Exact sign-flip enumeration exceeds `max_exact`; request a finite Monte Carlo sample.",
         call. = FALSE)
  }
  if (!exact_requested &&
      (length(permutations) != 1L || !is.finite(permutations) ||
       permutations < 0 || permutations != as.integer(permutations))) {
    stop("`permutations` must be a non-negative integer or \"all\".",
         call. = FALSE)
  }
  permutations <- if (exact_requested) NA_integer_ else as.integer(permutations)
  enumerate <- exact_requested ||
    (permutations > 0L && is.finite(attainable) &&
       permutations >= attainable - 1L && attainable <= max_exact)
  signs <- if (enumerate) {
    as.matrix(expand.grid(
      rep(list(c(-1, 1)), units), KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    ))
  } else if (permutations > 0L) {
    matrix(
      sample(c(-1, 1), permutations * units, replace = TRUE),
      nrow = permutations, ncol = units
    )
  } else matrix(numeric(), nrow = 0L, ncol = units)
  storage.mode(signs) <- "double"
  list(signs = signs, exact = enumerate, attainable = attainable)
}

.validate_random_slope_design <- function(design, tol = 1e-10) {
  if (!inherits(design, "permanova_design") || is.null(design$G)) {
    stop("`design` must be a compiled Design Object containing `distance`.",
         call. = FALSE)
  }
  if (length(design$random_terms) != 1L ||
      design$random_terms[[1L]]$origin != "slope") {
    stop("Random-slope inference requires exactly one `(1 + condition | subject)` term.",
         call. = FALSE)
  }
  specification <- design$random_terms[[1L]]
  slope_name <- specification$slope
  subject_name <- specification$group
  fixed_variables <- all.vars(design$fixed_formula)
  if (length(design$fixed_effects) != 1L ||
      length(fixed_variables) != 1L ||
      !identical(fixed_variables[[1L]], slope_name) ||
      !identical(.canonical_term(design$fixed_effects[[1L]]), slope_name)) {
    stop("Phase 27 requires the random-slope condition as the sole fixed main effect.",
         call. = FALSE)
  }
  raw_condition <- design$data[[slope_name]]
  if (is.numeric(raw_condition) && !is.factor(raw_condition)) {
    stop("The random-slope variable `", slope_name, "` is a continuous ",
         "covariate; use `permanova_random_slope_continuous()` (or the ",
         "`permanova()` front door, which routes it automatically).",
         call. = FALSE)
  }
  condition <- droplevels(as.factor(raw_condition))
  subject <- droplevels(factor(design$data[[subject_name]]))
  if (!is.factor(condition) || nlevels(condition) != 2L) {
    stop("Phase 27 requires one two-level categorical random-slope condition.",
         call. = FALSE)
  }
  if (nlevels(subject) < 3L) {
    stop("At least three complete subjects are required for random-slope inference.",
         call. = FALSE)
  }
  if (length(design$estimability$aliased) || any(!design$estimability$terms)) {
    stop("The fixed random-slope condition is rank deficient or aliased.",
         call. = FALSE)
  }
  cells <- interaction(subject, condition, drop = TRUE, lex.order = TRUE)
  incidence <- table(subject, condition)
  if (any(incidence < 1L)) {
    stop("Every subject must contain both random-slope condition levels.",
         call. = FALSE)
  }
  cell_sizes <- as.integer(table(cells))
  if (length(cell_sizes) != nlevels(subject) * 2L ||
      length(unique(cell_sizes)) != 1L || cell_sizes[[1L]] < 2L) {
    stop(
      paste(
        "Phase 27 requires equal replication of at least two observations",
        "in every subject-by-condition cell."
      ),
      call. = FALSE
    )
  }
  subjects <- levels(subject)
  conditions <- levels(condition)
  midpoint <- difference <- matrix(0, nrow = length(subjects), ncol = design$n)
  cell_blocks <- vector("list", length(subjects) * 2L)
  cell_index <- 0L
  for (subject_index in seq_along(subjects)) {
    for (condition_index in seq_along(conditions)) {
      observations <- which(
        subject == subjects[[subject_index]] &
          condition == conditions[[condition_index]]
      )
      cell_index <- cell_index + 1L
      cell_blocks[[cell_index]] <- observations
      midpoint[subject_index, observations] <-
        1 / (2 * length(observations))
      difference[subject_index, observations] <-
        c(-1, 1)[[condition_index]] / length(observations)
    }
  }
  K_midpoint <- midpoint %*% design$G %*% t(midpoint)
  K_difference <- difference %*% design$G %*% t(difference)
  K_cross <- midpoint %*% design$G %*% t(difference)
  K_midpoint <- (K_midpoint + t(K_midpoint)) / 2
  K_difference <- (K_difference + t(K_difference)) / 2
  number_subjects <- length(subjects)
  H_mean <- matrix(1 / number_subjects, number_subjects, number_subjects)
  H_subject_residual <- diag(number_subjects) - H_mean
  cell_factor <- factor(cells)
  H_cell <- projection_matrix(
    stats::model.matrix(~ 0 + cell_factor), tol = tol
  )$H
  within_df <- design$n - nlevels(cell_factor)
  if (within_df <= 0L) {
    stop("No within-cell residual degrees of freedom remain.", call. = FALSE)
  }
  list(
    slope_name = slope_name, subject_name = subject_name,
    condition = condition, subject = subject,
    conditions = conditions, subjects = subjects,
    replicates = cell_sizes[[1L]], cell_blocks = cell_blocks,
    midpoint_membership = midpoint, difference_membership = difference,
    K_midpoint = K_midpoint, K_difference = K_difference,
    K_cross = K_cross, H_mean = H_mean,
    H_subject_residual = H_subject_residual,
    numerator_df = 1L, denominator_df = number_subjects - 1L,
    H_cell = H_cell, H_within = diag(design$n) - H_cell,
    within_df = within_df
  )
}

#' Estimate a distance-based random intercept--slope covariance
#'
#' Returns trace covariance components for the subject midpoint (random
#' intercept at centred condition coding) and the second-minus-first condition
#' difference (random slope). Equal-replication within-cell noise is removed by
#' method of moments. Negative unrestricted components are retained; optional
#' eigenvalue truncation is a separate sensitivity.
#'
#' @param design Canonical Phase 27 random-slope Design Object.
#' @param truncate Also return a positive-semidefinite eigenvalue-truncated
#'   sensitivity matrix.
#' @param tol Numerical tolerance.
#' @return Random-coefficient covariance estimates and diagnostics.
#' @export
estimate_random_slope_covariance <- function(
    design, truncate = FALSE, tol = 1e-10) {
  structure <- .validate_random_slope_design(design, tol = tol)
  center <- structure$H_subject_residual
  subjects_df <- structure$denominator_df
  residual_ss <- trace_hgh(structure$H_within, design$G)
  residual_ms <- residual_ss / structure$within_df
  midpoint_variance <- trace_product(
    center, structure$K_midpoint
  ) / subjects_df
  slope_variance <- trace_product(
    center, structure$K_difference
  ) / subjects_df
  intercept_slope_covariance <- trace_product(
    center, structure$K_cross
  ) / subjects_df
  unrestricted <- matrix(c(
    midpoint_variance - residual_ms / (2 * structure$replicates),
    intercept_slope_covariance,
    intercept_slope_covariance,
    slope_variance - 2 * residual_ms / structure$replicates
  ), 2L, 2L, byrow = TRUE, dimnames = list(
    c("random_intercept", "random_slope"),
    c("random_intercept", "random_slope")
  ))
  truncated <- NULL
  if (truncate) {
    decomposition <- eigen(unrestricted, symmetric = TRUE)
    truncated <- decomposition$vectors %*%
      diag(pmax(decomposition$values, 0), 2L) %*%
      t(decomposition$vectors)
    dimnames(truncated) <- dimnames(unrestricted)
  }
  correlation <- if (all(diag(unrestricted) > tol)) {
    unrestricted[1L, 2L] / sqrt(prod(diag(unrestricted)))
  } else NA_real_
  list(
    unrestricted = unrestricted, truncated = truncated,
    residual_trace_variance = residual_ms,
    empirical_midpoint_variance = midpoint_variance,
    empirical_difference_variance = slope_variance,
    intercept_slope_correlation = correlation,
    parameterization = paste0(
      "random intercept is the midpoint of `", structure$conditions[[1L]],
      "` and `", structure$conditions[[2L]],
      "`; random slope is `", structure$conditions[[2L]], " - ",
      structure$conditions[[1L]], "`"
    ),
    warnings = c(
      if (any(diag(unrestricted) < -tol)) {
        "Negative unrestricted random-coefficient variance components were retained."
      } else character(),
      if (is.finite(correlation) && abs(correlation) > 1 + sqrt(tol)) {
        "The unrestricted trace correlation lies outside [-1, 1]."
      } else character()
    ),
    method = "balanced two-condition trace method of moments"
  )
}

#' Balanced two-condition random-slope PERMANOVA
#'
#' Tests the population mean within-subject condition difference for a
#' correlated `(1 + condition | subject)` model. Complete difference vectors
#' are sign-flipped by subject. Exactness is conditional on central symmetry of
#' subject difference vectors under the fixed-effect null; this is stated in
#' every result and is not presented as unrestricted label exchangeability.
#'
#' @param design Canonical Phase 27 random-slope Design Object.
#' @param permutations A non-negative integer or `"all"`.
#' @param seed Optional sign-flip seed.
#' @param keep_permutations Retain sign-flip statistics.
#' @param truncate Also return a positive-semidefinite covariance sensitivity.
#' @param max_exact Maximum exact sign-flip space.
#' @param tol Numerical tolerance.
#' @return An object of class `permanova_random_slope`.
#' @export
permanova_random_slope <- function(
    design, permutations = 999L, seed = NULL,
    keep_permutations = FALSE, truncate = FALSE,
    max_exact = 100000L, tol = 1e-10) {
  structure <- .validate_random_slope_design(design, tol = tol)
  numerator_ss <- trace_hgh(structure$H_mean, structure$K_difference)
  denominator_ss <- trace_hgh(
    structure$H_subject_residual, structure$K_difference
  )
  denominator_ms <- denominator_ss / structure$denominator_df
  if (!is.finite(denominator_ms) || abs(denominator_ms) <= tol) {
    stop("The between-subject random-slope denominator is numerically zero.",
         call. = FALSE)
  }
  observed <- numerator_ss / denominator_ms
  if (!is.null(seed)) set.seed(seed)
  plan <- .sign_flip_plan(
    length(structure$subjects), permutations, max_exact
  )
  statistics <- vapply(seq_len(nrow(plan$signs)), function(index) {
    signs <- plan$signs[index, ]
    permuted <- structure$K_difference * tcrossprod(signs)
    denominator <- trace_hgh(
      structure$H_subject_residual, permuted
    ) / structure$denominator_df
    if (!is.finite(denominator) || abs(denominator) <= tol) {
      stop("A sign-flipped random-slope denominator is numerically zero.",
           call. = FALSE)
    }
    trace_hgh(structure$H_mean, permuted) / denominator
  }, numeric(1))
  p_value <- NA_real_
  if (length(statistics)) {
    extreme <- sum(statistics >= observed - tol, na.rm = TRUE)
    p_value <- if (plan$exact) {
      extreme / length(statistics)
    } else (extreme + 1) / (length(statistics) + 1)
  }
  parametric_p <- if (isTRUE(design$spectral$euclidean)) {
    stats::pf(observed, 1L, structure$denominator_df, lower.tail = FALSE)
  } else NA_real_
  notes <- paste(
    "Subject sign flips require central symmetry of complete condition-difference",
    "vectors under the population fixed-effect null."
  )
  if (is.finite(plan$attainable) && plan$attainable < 100) {
    notes <- c(notes, paste0(
      "Only ", plan$attainable, " subject sign assignments are attainable."
    ))
  }
  if (!isTRUE(design$spectral$euclidean)) {
    notes <- c(
      notes,
      "Empirical inference retained the complete non-Euclidean Gower spectrum; the parametric p-value was suppressed."
    )
  }
  test <- list(
    term = structure$slope_name,
    hypothesis = paste0(
      "the population mean subject difference `",
      structure$conditions[[2L]], " - ", structure$conditions[[1L]],
      "` is zero"
    ),
    estimand = "equally weighted complete-subject condition-difference vectors",
    numerator = list(
      ss = numerator_ss, df = 1L, ms = numerator_ss,
      projector = structure$H_mean
    ),
    denominator = list(
      expression = "MS[between-subject random-slope differences]",
      ss = denominator_ss, df = structure$denominator_df,
      ms = denominator_ms, projector = structure$H_subject_residual
    ),
    statistic = observed, p_value = p_value,
    p_value_parametric = parametric_p,
    exchangeable_units = structure$subject_name,
    permutation_strata = "independent whole-subject difference vectors",
    permutation_algorithm = "whole-subject sign flips of complete condition-difference vectors",
    exchangeability_assumption = "central symmetry of subject difference vectors under the fixed-effect null",
    attainable_permutations = plan$attainable,
    permutations_used = nrow(plan$signs),
    exact_permutation = plan$exact && nrow(plan$signs) > 0L,
    approximation_status = if (!nrow(plan$signs)) {
      "not tested"
    } else if (plan$exact) {
      "exhaustive whole-subject sign-flip inference"
    } else "Monte Carlo whole-subject sign-flip inference",
    warnings = notes,
    permutation_statistics = if (keep_permutations) statistics else NULL
  )
  class(test) <- "permanova_test_spec"
  table <- data.frame(
    term = c(structure$slope_name, "Between-subject random slope"),
    df = c(1L, structure$denominator_df),
    ss = c(numerator_ss, denominator_ss),
    ms = c(numerator_ss, denominator_ms),
    pseudo_f = c(observed, NA_real_),
    p_sign = c(p_value, NA_real_),
    p_parametric = c(parametric_p, NA_real_)
  )
  out <- list(
    call = match.call(), design = design, table = table, test = test,
    covariance = estimate_random_slope_covariance(
      design, truncate = truncate, tol = tol
    ),
    structure = list(
      condition = structure$slope_name, subject = structure$subject_name,
      condition_levels = structure$conditions,
      subjects = length(structure$subjects),
      replicates_per_cell = structure$replicates
    ),
    estimand = test$estimand,
    supported_design = paste(
      "balanced complete two-condition repeated measures with one correlated",
      "random intercept and slope; central-symmetry sign-flip inference"
    ),
    permutation = list(
      algorithm = test$permutation_algorithm,
      assumption = test$exchangeability_assumption,
      exact = test$exact_permutation, attainable = plan$attainable,
      used = nrow(plan$signs), seed = seed
    )
  )
  class(out) <- "permanova_random_slope"
  out
}

#' @export
print.permanova_random_slope <- function(x, ...) {
  cat("Balanced two-condition random-slope PERMANOVA\n")
  cat("Estimand:", x$estimand, "\n")
  cat("Sign-flip assumption:", x$permutation$assumption, "\n")
  print(x$table, row.names = FALSE)
  cat("Unrestricted random-coefficient trace covariance:\n")
  print(x$covariance$unrestricted)
  invisible(x)
}
