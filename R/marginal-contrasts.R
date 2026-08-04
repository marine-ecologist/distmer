.validate_marginal_design <- function(design, factor, tol = 1e-10) {
  if (!inherits(design, "permanova_design") || is.null(design$G)) {
    stop("`design` must be a compiled Design Object containing `distance`.",
         call. = FALSE)
  }
  if (length(design$random_terms)) {
    stop("Marginal contrasts support fixed-effect models only.", call. = FALSE)
  }
  variables <- all.vars(design$fixed_formula)
  if (!length(variables) || any(!vapply(design$data[variables], is.factor,
                                        logical(1)))) {
    stop("Marginal contrasts currently require every fixed predictor to be categorical.",
         call. = FALSE)
  }
  if (length(factor) != 1L || !is.character(factor) ||
      !factor %in% variables) {
    stop("`factor` must name one categorical predictor in the fixed model.",
         call. = FALSE)
  }
  if (!factor %in% design$fixed_effects) {
    stop("The target factor must occur as an explicit fixed main effect.",
         call. = FALSE)
  }
  if (length(design$estimability$aliased) ||
      design$estimability$fixed_rank != ncol(design$X)) {
    stop("Marginal contrasts require a full-rank fixed model; aliased reference-grid predictions are unsupported.",
         call. = FALSE)
  }
  denominator_df <- design$n - design$estimability$fixed_rank
  if (denominator_df <= 0L) {
    stop("No full-model residual degrees of freedom remain.", call. = FALSE)
  }
  list(
    variables = variables, target = factor,
    target_levels = levels(droplevels(design$data[[factor]])),
    nuisance = setdiff(variables, factor),
    H_full = design$projections$fixed$H,
    H_denominator = diag(design$n) - design$projections$fixed$H,
    denominator_df = denominator_df
  )
}

.factor_grid <- function(data, variables) {
  if (!length(variables)) {
    return(data.frame(.row = 1L)[, FALSE, drop = FALSE])
  }
  grid <- expand.grid(
    lapply(data[variables], function(x) levels(droplevels(x))),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  names(grid) <- variables
  for (variable in variables) {
    grid[[variable]] <- factor(
      grid[[variable]], levels = levels(droplevels(data[[variable]]))
    )
  }
  grid
}

.observed_nuisance_grid <- function(data, variables) {
  if (!length(variables)) {
    return(list(
      grid = data.frame(.row = 1L)[, FALSE, drop = FALSE], weights = 1
    ))
  }
  key <- interaction(data[variables], drop = TRUE, lex.order = TRUE)
  first <- !duplicated(key)
  grid <- data[first, variables, drop = FALSE]
  weights <- as.numeric(table(key))[match(as.character(key[first]), names(table(key)))]
  weights <- weights / sum(weights)
  ordering <- do.call(order, lapply(grid, unclass))
  grid <- grid[ordering, , drop = FALSE]
  weights <- weights[ordering]
  rownames(grid) <- NULL
  list(grid = grid, weights = weights)
}

.restore_factor_contract <- function(grid, design, variables) {
  for (variable in variables) {
    original <- design$data[[variable]]
    grid[[variable]] <- factor(grid[[variable]], levels = levels(original))
    original_contrasts <- stats::contrasts(original)
    if (!is.null(original_contrasts)) {
      stats::contrasts(grid[[variable]]) <- original_contrasts
    }
  }
  grid
}

#' Derive fixed-model marginal means
#'
#' Constructs coefficient-space rows for every level of a target factor in a
#' full-rank categorical fixed model. Equal weighting averages the complete
#' nuisance-factor reference grid. Proportional weighting averages observed
#' nuisance combinations by their joint frequency.
#'
#' @param design Canonical fixed-effect Design Object.
#' @param factor Target categorical predictor.
#' @param weights `"equal"` or `"proportional"` nuisance-cell weights.
#' @param tol Numerical tolerance.
#' @return Marginal-mean coefficient rows and reference-grid metadata.
#' @export
derive_marginal_means <- function(
    design, factor, weights = c("equal", "proportional"), tol = 1e-10) {
  weights <- match.arg(weights)
  structure <- .validate_marginal_design(design, factor, tol = tol)
  nuisance <- if (weights == "equal") {
    list(
      grid = .factor_grid(design$data, structure$nuisance),
      weights = NULL
    )
  } else {
    .observed_nuisance_grid(design$data, structure$nuisance)
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
  for (index in seq_along(structure$target_levels)) {
    grid <- nuisance$grid
    grid[[factor]] <- factor(
      rep(structure$target_levels[[index]], nrow(grid)),
      levels = levels(design$data[[factor]])
    )
    grid <- grid[, structure$variables, drop = FALSE]
    grid <- .restore_factor_contract(grid, design, structure$variables)
    frame <- stats::model.frame(
      stats::terms(design$fixed_formula), data = grid,
      na.action = stats::na.fail, drop.unused.levels = FALSE
    )
    matrix <- stats::model.matrix(
      stats::terms(design$fixed_formula), frame
    )
    if (!identical(colnames(matrix), colnames(design$X))) {
      stop("Reference-grid model columns do not match the compiled design.",
           call. = FALSE)
    }
    coefficient_rows[index, ] <- colSums(matrix * nuisance$weights)
    reference_grids[[index]] <- grid
  }
  list(
    factor = factor, levels = structure$target_levels,
    weighting = weights, coefficient_rows = coefficient_rows,
    nuisance_factors = structure$nuisance,
    nuisance_grid = nuisance$grid, nuisance_weights = nuisance$weights,
    reference_grids = reference_grids,
    estimable = TRUE,
    definition = if (weights == "equal") {
      "equal weight over the complete nuisance-factor reference grid"
    } else {
      "weight over observed joint nuisance combinations in proportion to frequency"
    }
  )
}

.linear_contrast_projector <- function(X, coefficient_row, H_full, tol) {
  information <- crossprod(X)
  direction <- as.vector(X %*% qr.solve(
    information, as.numeric(coefficient_row), tol = tol
  ))
  length_squared <- sum(direction ^ 2)
  if (!is.finite(length_squared) || length_squared <= tol) {
    stop("The requested marginal contrast has numerically zero estimable length.",
         call. = FALSE)
  }
  H <- tcrossprod(direction) / length_squared
  H_reduced <- H_full - H
  error <- max(
    abs(H %*% H - H),
    abs(H_reduced %*% H_reduced - H_reduced),
    abs(H %*% H_reduced)
  )
  if (error > max(100 * tol, sqrt(.Machine$double.eps))) {
    stop("The marginal contrast does not form an estimable one-df projector.",
         call. = FALSE)
  }
  list(H = H, H_reduced = H_reduced, direction = direction, error = error)
}

#' Multifactor marginal contrasts
#'
#' Tests contrasts among explicit marginal means from a full-rank categorical
#' fixed model. Each contrast has its own nuisance-retaining reduced model and
#' unrestricted Freedman--Lane residual permutation.
#'
#' @inheritParams derive_marginal_means
#' @param contrasts Numeric level contrast vector, row matrix, or named list.
#' @param adjust Multiplicity adjustment.
#' @param permutations A non-negative integer or `"all"`.
#' @param seed Optional random seed.
#' @param keep_permutations Retain permutation statistics.
#' @param max_exact Maximum exact permutation space.
#' @return An object of class `permanova_marginal_contrasts`.
#' @export
permanova_marginal_contrasts <- function(
    design, factor, contrasts, weights = c("equal", "proportional"),
    adjust = c("holm", "bonferroni", "BH", "BY", "none"),
    permutations = 999L, seed = NULL, keep_permutations = FALSE,
    max_exact = 100000L, tol = 1e-10) {
  weights <- match.arg(weights)
  adjust <- match.arg(adjust)
  structure <- .validate_marginal_design(design, factor, tol = tol)
  means <- derive_marginal_means(design, factor, weights, tol)
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
      structure$H_denominator, structure$denominator_df, plan, tol,
      label = paste0("Permuting marginal contrast `",
                     rownames(contrast_matrix)[[index]], "`")
    )
    test <- list(
      term = factor, contrast = rownames(contrast_matrix)[[index]],
      coefficients = coefficients,
      marginal_coefficient_row = setNames(coefficient_row, colnames(design$X)),
      weighting = weights, weighting_definition = means$definition,
      hypothesis = paste0("the `", weights, "` marginal contrast `",
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
      permutation_algorithm = "marginal-contrast-specific Freedman-Lane residual permutation",
      attainable_permutations = plan$attainable,
      permutations_used = nrow(plan$orders),
      exact_permutation = plan$exact && nrow(plan$orders) > 0L,
      approximation_status = if (!nrow(plan$orders)) "not tested" else if (plan$exact) {
        "exhaustive marginal-contrast residual permutation"
      } else "Monte Carlo marginal-contrast residual permutation",
      projector_error = projector$error,
      warnings = if (!isTRUE(design$spectral$euclidean)) {
        "Empirical inference retained the complete non-Euclidean Gower spectrum."
      } else character(),
      permutation_statistics = if (keep_permutations) result$statistics else NULL
    )
    class(test) <- "permanova_test_spec"
    test
  })
  names(tests) <- rownames(contrast_matrix)
  raw_p <- vapply(tests, `[[`, numeric(1), "p_value")
  adjusted <- if (all(is.na(raw_p))) raw_p else stats::p.adjust(raw_p, adjust)
  for (index in seq_along(tests)) tests[[index]]$p_adjusted <- adjusted[[index]]
  table <- data.frame(
    contrast = names(tests), weighting = weights,
    df = 1L,
    ss = vapply(tests, function(x) x$numerator$ss, numeric(1)),
    denominator_df = structure$denominator_df,
    pseudo_f = vapply(tests, `[[`, numeric(1), "statistic"),
    p_perm = raw_p, p_adjusted = adjusted, adjustment = adjust,
    stringsAsFactors = FALSE
  )
  out <- list(
    call = match.call(), design = design, factor = factor,
    marginal_means = means, contrast_matrix = contrast_matrix,
    table = table, tests = tests, adjustment = adjust,
    supported_design = paste(
      "full-rank categorical fixed model with explicit", weights,
      "marginal-mean weighting"
    ),
    permutation = list(
      algorithm = "contrast-specific unrestricted Freedman-Lane residual permutation",
      exact = plan$exact && nrow(plan$orders) > 0L,
      attainable = plan$attainable, used = nrow(plan$orders), seed = seed
    )
  )
  class(out) <- "permanova_marginal_contrasts"
  out
}

#' Pairwise multifactor marginal comparisons
#'
#' @inheritParams permanova_marginal_contrasts
#' @return An object of class `permanova_marginal_contrasts`.
#' @export
permanova_marginal_pairwise <- function(
    design, factor, weights = c("equal", "proportional"),
    adjust = c("holm", "bonferroni", "BH", "BY", "none"),
    permutations = 999L, seed = NULL, keep_permutations = FALSE,
    max_exact = 100000L, tol = 1e-10) {
  structure <- .validate_marginal_design(design, factor, tol)
  pairs <- utils::combn(structure$target_levels, 2L, simplify = FALSE)
  contrast_matrix <- matrix(
    0, nrow = length(pairs), ncol = length(structure$target_levels),
    dimnames = list(NULL, structure$target_levels)
  )
  for (index in seq_along(pairs)) {
    contrast_matrix[index, pairs[[index]][[1L]]] <- 1
    contrast_matrix[index, pairs[[index]][[2L]]] <- -1
  }
  rownames(contrast_matrix) <- vapply(
    pairs, function(x) paste(x, collapse = " - "), character(1)
  )
  out <- permanova_marginal_contrasts(
    design, factor, contrast_matrix, weights = match.arg(weights),
    adjust = match.arg(adjust), permutations = permutations, seed = seed,
    keep_permutations = keep_permutations, max_exact = max_exact, tol = tol
  )
  out$call <- match.call()
  out
}

#' @export
print.permanova_marginal_contrasts <- function(x, ...) {
  cat("Multifactor marginal distance contrasts\n")
  cat("Target:", x$factor, " Weighting:", x$marginal_means$weighting, "\n")
  print(x$table, row.names = FALSE)
  invisible(x)
}
