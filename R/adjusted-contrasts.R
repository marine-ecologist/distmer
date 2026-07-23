.normalise_reference_values <- function(at, continuous) {
  if (!length(continuous)) {
    if (!is.null(at) && length(at)) {
      stop("`at` must be empty when the fixed model has no continuous nuisance variables.",
           call. = FALSE)
    }
    return(setNames(numeric(), character()))
  }
  if (is.list(at) && !is.data.frame(at)) at <- unlist(at, use.names = TRUE)
  if (!is.numeric(at) || is.null(names(at)) || any(!nzchar(names(at))) ||
      anyDuplicated(names(at)) || !setequal(names(at), continuous)) {
    stop(
      "`at` must be a named numeric scalar for every continuous nuisance variable and no others: ",
      paste(continuous, collapse = ", "), ".", call. = FALSE
    )
  }
  at <- at[continuous]
  if (length(at) != length(continuous) || any(!is.finite(at))) {
    stop("Every continuous reference value must be one finite numeric scalar.",
         call. = FALSE)
  }
  at
}

.validate_adjusted_design <- function(design, factor, at, tol = 1e-10) {
  if (length(tol) != 1L || !is.numeric(tol) || !is.finite(tol) || tol <= 0) {
    stop("`tol` must be one finite positive number.", call. = FALSE)
  }
  if (!inherits(design, "permanova_design") || is.null(design$G)) {
    stop("`design` must be a compiled Design Object containing `distance`.",
         call. = FALSE)
  }
  if (length(design$random_terms)) {
    stop("Covariate-adjusted contrasts support fixed-effect models only.",
         call. = FALSE)
  }
  variables <- all.vars(design$fixed_formula)
  if (length(factor) != 1L || !is.character(factor) ||
      !factor %in% variables || !is.factor(design$data[[factor]])) {
    stop("`factor` must name one categorical predictor in the fixed model.",
         call. = FALSE)
  }
  if (!factor %in% design$fixed_effects) {
    stop("The target factor must occur as an explicit fixed main effect.",
         call. = FALSE)
  }
  supported <- vapply(design$data[variables], function(value) {
    is.factor(value) || (is.numeric(value) && !is.object(value))
  }, logical(1))
  if (any(!supported)) {
    stop("Fixed predictors must be factors or ordinary numeric covariates.",
         call. = FALSE)
  }
  continuous <- setdiff(
    variables[vapply(design$data[variables], is.numeric, logical(1))], factor
  )
  reference_values <- .normalise_reference_values(at, continuous)
  if (length(design$estimability$aliased) ||
      design$estimability$fixed_rank != ncol(design$X)) {
    stop("Covariate-adjusted contrasts require a full-rank fixed model.",
         call. = FALSE)
  }
  singular_values <- svd(design$X, nu = 0L, nv = 0L)$d
  condition_ratio <- min(singular_values) / max(singular_values)
  if (!is.finite(condition_ratio) || condition_ratio <= sqrt(tol)) {
    stop(
      paste0(
        "Covariate-adjusted contrasts require a numerically full-rank fixed ",
        "model; the fixed-model matrix is ill-conditioned at the requested ",
        "tolerance."
      ),
      call. = FALSE
    )
  }
  denominator_df <- design$n - design$estimability$fixed_rank
  if (denominator_df <= 0L) {
    stop("No full-model residual degrees of freedom remain.", call. = FALSE)
  }
  list(
    variables = variables, target = factor,
    target_levels = levels(droplevels(design$data[[factor]])),
    continuous = continuous,
    categorical = setdiff(
      variables[vapply(design$data[variables], is.factor, logical(1))],
      factor
    ),
    reference_values = reference_values,
    fixed_model_condition_ratio = condition_ratio,
    H_full = design$projections$fixed$H,
    H_denominator = diag(design$n) - design$projections$fixed$H,
    denominator_df = denominator_df
  )
}

#' Derive explicitly covariate-adjusted means
#'
#' Constructs coefficient rows for a target factor while holding every
#' continuous nuisance variable at a user-supplied value. Categorical nuisance
#' factors are averaged either over their complete reference grid or over
#' observed joint frequencies. No covariate reference value is inferred.
#'
#' @param design Canonical fixed-effect Design Object.
#' @param factor Target categorical predictor.
#' @param at Named numeric vector or list giving one finite scalar for every
#'   continuous nuisance variable.
#' @param weights `"equal"` or `"proportional"` categorical nuisance weights.
#' @param tol Numerical tolerance.
#' @return Adjusted-mean coefficient rows and reference-grid metadata.
#' @export
derive_adjusted_means <- function(
    design, factor, at = NULL,
    weights = c("equal", "proportional"), tol = 1e-10) {
  weights <- match.arg(weights)
  structure <- .validate_adjusted_design(design, factor, at, tol)
  nuisance <- if (weights == "equal") {
    list(grid = .factor_grid(design$data, structure$categorical),
         weights = NULL)
  } else {
    .observed_nuisance_grid(design$data, structure$categorical)
  }
  if (is.null(nuisance$weights)) {
    nuisance$weights <- rep(1 / nrow(nuisance$grid), nrow(nuisance$grid))
  }
  coefficient_rows <- matrix(
    0, nrow = length(structure$target_levels), ncol = ncol(design$X),
    dimnames = list(structure$target_levels, colnames(design$X))
  )
  reference_grids <- vector("list", length(structure$target_levels))
  names(reference_grids) <- structure$target_levels
  factor_variables <- c(factor, structure$categorical)
  for (index in seq_along(structure$target_levels)) {
    grid <- nuisance$grid
    for (variable in structure$continuous) {
      grid[[variable]] <- rep(structure$reference_values[[variable]], nrow(grid))
    }
    grid[[factor]] <- factor(
      rep(structure$target_levels[[index]], nrow(grid)),
      levels = levels(design$data[[factor]])
    )
    grid <- grid[, structure$variables, drop = FALSE]
    grid <- .restore_factor_contract(
      grid, design, intersect(structure$variables, factor_variables)
    )
    frame <- stats::model.frame(
      stats::terms(design$fixed_formula), data = grid,
      na.action = stats::na.fail, drop.unused.levels = FALSE
    )
    matrix <- stats::model.matrix(stats::terms(design$fixed_formula), frame)
    if (!identical(colnames(matrix), colnames(design$X))) {
      stop("Adjusted reference-grid columns do not match the compiled design.",
           call. = FALSE)
    }
    coefficient_rows[index, ] <- colSums(matrix * nuisance$weights)
    reference_grids[[index]] <- grid
  }
  list(
    factor = factor, levels = structure$target_levels,
    weighting = weights, coefficient_rows = coefficient_rows,
    categorical_nuisance = structure$categorical,
    continuous_nuisance = structure$continuous,
    reference_values = structure$reference_values,
    fixed_model_condition_ratio = structure$fixed_model_condition_ratio,
    nuisance_grid = nuisance$grid, nuisance_weights = nuisance$weights,
    reference_grids = reference_grids, estimable = TRUE,
    definition = paste0(
      if (weights == "equal") {
        "equal weight over the categorical nuisance reference grid"
      } else {
        "observed-frequency weight over joint categorical nuisance cells"
      },
      if (length(structure$continuous)) paste0(
        "; continuous covariates fixed explicitly at ",
        paste(names(structure$reference_values),
              format(structure$reference_values), sep = "=", collapse = ", ")
      ) else ""
    )
  )
}

#' Explicitly covariate-adjusted distance contrasts
#'
#' Tests contrasts among adjusted means with contrast-specific Freedman--Lane
#' residual permutations and the common full-model residual denominator.
#'
#' @inheritParams derive_adjusted_means
#' @param contrasts Numeric level contrast vector, row matrix, or named list.
#' @param adjust Multiplicity adjustment.
#' @param permutations A non-negative integer or `"all"`.
#' @param seed Optional random seed.
#' @param keep_permutations Retain permutation statistics.
#' @param max_exact Maximum exact permutation space.
#' @return An object of class `permanova_adjusted_contrasts`.
#' @export
permanova_adjusted_contrasts <- function(
    design, factor, contrasts, at = NULL,
    weights = c("equal", "proportional"),
    adjust = c("holm", "bonferroni", "BH", "BY", "none"),
    permutations = 999L, seed = NULL, keep_permutations = FALSE,
    max_exact = 100000L, tol = 1e-10) {
  weights <- match.arg(weights)
  adjust <- match.arg(adjust)
  structure <- .validate_adjusted_design(design, factor, at, tol)
  means <- derive_adjusted_means(design, factor, at, weights, tol)
  contrast_matrix <- .as_contrast_matrix(
    contrasts, structure$target_levels, tol
  )
  if (!is.null(seed)) set.seed(seed)
  plan <- .permutation_orders(design$n, permutations, max_exact)
  tests <- lapply(seq_len(nrow(contrast_matrix)), function(index) {
    coefficients <- contrast_matrix[index, ]
    coefficient_row <- as.numeric(coefficients %*% means$coefficient_rows)
    projector <- .linear_contrast_projector(
      design$X, coefficient_row, structure$H_full, tol
    )
    result <- .unit_freedman_lane(
      design$G, projector$H, 1L, projector$H_reduced,
      structure$H_denominator, structure$denominator_df, plan, tol
    )
    notes <- if (!isTRUE(design$spectral$euclidean)) {
      "Empirical inference retained the complete non-Euclidean Gower spectrum."
    } else character()
    test <- list(
      term = factor, contrast = rownames(contrast_matrix)[[index]],
      coefficients = coefficients,
      adjusted_coefficient_row = setNames(coefficient_row, colnames(design$X)),
      reference_values = structure$reference_values,
      weighting = weights, weighting_definition = means$definition,
      hypothesis = paste0("the explicitly adjusted contrast `",
                          rownames(contrast_matrix)[[index]], "` is zero"),
      numerator = list(
        ss = result$numerator_ss, df = 1L, ms = result$numerator_ms,
        projector = projector$H
      ),
      denominator = list(
        expression = "MS[full-model Residual]",
        ss = result$denominator_ss, df = structure$denominator_df,
        ms = result$denominator_ms,
        projector = structure$H_denominator
      ),
      statistic = result$observed, p_value = result$p_value,
      p_adjusted = NA_real_, adjustment = adjust,
      exchangeable_units = "observations", permutation_strata = "unrestricted",
      permutation_algorithm = "adjusted-contrast-specific Freedman-Lane residual permutation",
      attainable_permutations = plan$attainable,
      permutations_used = nrow(plan$orders),
      exact_permutation = plan$exact && nrow(plan$orders) > 0L,
      approximation_status = if (!nrow(plan$orders)) {
        "not tested"
      } else if (plan$exact) {
        "exhaustive adjusted-contrast residual permutation"
      } else "Monte Carlo adjusted-contrast residual permutation",
      warnings = notes,
      permutation_statistics = if (keep_permutations) result$statistics else NULL
    )
    class(test) <- "permanova_test_spec"
    test
  })
  names(tests) <- rownames(contrast_matrix)
  raw <- vapply(tests, `[[`, numeric(1), "p_value")
  adjusted <- if (all(is.na(raw))) raw else stats::p.adjust(raw, adjust)
  for (index in seq_along(tests)) tests[[index]]$p_adjusted <- adjusted[[index]]
  table <- data.frame(
    contrast = names(tests), weighting = weights, df = 1L,
    ss = vapply(tests, function(test) test$numerator$ss, numeric(1)),
    denominator_df = structure$denominator_df,
    pseudo_f = vapply(tests, `[[`, numeric(1), "statistic"),
    p_perm = raw, p_adjusted = adjusted, adjustment = adjust,
    stringsAsFactors = FALSE
  )
  out <- list(
    call = match.call(), design = design, factor = factor,
    levels = structure$target_levels, contrast_matrix = contrast_matrix,
    adjusted_means = means, reference_values = structure$reference_values,
    table = table, tests = tests, adjustment = adjust,
    family = "user-specified explicitly covariate-adjusted contrasts",
    supported_design = "full-rank fixed model with categorical factors and explicitly referenced numeric covariates",
    permutation = list(
      algorithm = "contrast-specific unrestricted Freedman-Lane residual permutation",
      exact = plan$exact && nrow(plan$orders) > 0L,
      attainable = plan$attainable, used = nrow(plan$orders), seed = seed
    )
  )
  class(out) <- "permanova_adjusted_contrasts"
  out
}

#' Pairwise explicitly covariate-adjusted distance contrasts
#'
#' @inheritParams permanova_adjusted_contrasts
#' @return An object of class `permanova_adjusted_contrasts`.
#' @export
permanova_adjusted_pairwise <- function(
    design, factor, at = NULL,
    weights = c("equal", "proportional"),
    adjust = c("holm", "bonferroni", "BH", "BY", "none"),
    permutations = 999L, seed = NULL, keep_permutations = FALSE,
    max_exact = 100000L, tol = 1e-10) {
  structure <- .validate_adjusted_design(design, factor, at, tol)
  pairs <- utils::combn(structure$target_levels, 2L, simplify = FALSE)
  contrasts <- matrix(
    0, nrow = length(pairs), ncol = length(structure$target_levels),
    dimnames = list(NULL, structure$target_levels)
  )
  for (index in seq_along(pairs)) {
    contrasts[index, pairs[[index]][[1L]]] <- 1
    contrasts[index, pairs[[index]][[2L]]] <- -1
  }
  rownames(contrasts) <- vapply(
    pairs, function(pair) paste(pair, collapse = " - "), character(1)
  )
  out <- permanova_adjusted_contrasts(
    design, factor, contrasts, at, weights = match.arg(weights),
    adjust = match.arg(adjust), permutations = permutations, seed = seed,
    keep_permutations = keep_permutations, max_exact = max_exact, tol = tol
  )
  out$call <- match.call()
  out$family <- "all pairwise explicitly covariate-adjusted contrasts"
  out
}

#' @export
print.permanova_adjusted_contrasts <- function(x, ...) {
  cat("Explicitly covariate-adjusted distance contrasts\n")
  cat("Reference values:",
      paste(names(x$reference_values), format(x$reference_values),
            sep = "=", collapse = ", "), "\n")
  cat("Family:", x$family, "\n")
  print(x$table, row.names = FALSE)
  invisible(x)
}
