.validate_contrast_design <- function(design, tol = 1e-10) {
  if (!inherits(design, "permanova_design") || is.null(design$G)) {
    stop("`design` must be a compiled Design Object containing `distance`.",
         call. = FALSE)
  }
  if (length(design$random_terms)) {
    stop(
      "Phase 13 contrasts support fixed one-factor designs only; mixed-model contrast denominators require a separate design-specific engine.",
      call. = FALSE
    )
  }
  fixed <- .one_categorical_fixed(design, "Phase 13 contrasts")
  group <- droplevels(fixed$value)
  counts <- table(group)
  if (any(counts < 1L)) {
    stop("Every factor level must contain observations.", call. = FALSE)
  }
  full <- design$projections$fixed
  if (full$rank != nlevels(group) || nrow(design$X) - full$rank <= 0L ||
      length(design$estimability$aliased)) {
    stop("The one-factor model must be full rank with positive within-group residual df.",
         call. = FALSE)
  }
  list(
    fixed = fixed, group = group, levels = levels(group), counts = counts,
    H_full = full$H, denominator_df = design$n - full$rank,
    H_denominator = diag(design$n) - full$H
  )
}

.as_contrast_matrix <- function(contrasts, levels, tol) {
  if (is.list(contrasts) && !is.data.frame(contrasts)) {
    if (!length(contrasts) || any(!vapply(contrasts, is.numeric, logical(1)))) {
      stop("`contrasts` must contain numeric contrast coefficients.",
           call. = FALSE)
    }
    matrix_value <- do.call(rbind, contrasts)
    if (!is.null(names(contrasts))) rownames(matrix_value) <- names(contrasts)
  } else if (is.numeric(contrasts) && is.null(dim(contrasts))) {
    matrix_value <- matrix(contrasts, nrow = 1L)
    if (!is.null(names(contrasts))) colnames(matrix_value) <- names(contrasts)
  } else if (is.matrix(contrasts) && is.numeric(contrasts)) {
    matrix_value <- contrasts
  } else {
    stop("`contrasts` must be a numeric vector, matrix, or named list.",
         call. = FALSE)
  }
  if (ncol(matrix_value) != length(levels)) {
    stop("Every contrast must supply one coefficient for each factor level.",
         call. = FALSE)
  }
  if (!is.null(colnames(matrix_value))) {
    if (!setequal(colnames(matrix_value), levels)) {
      stop("Contrast column names must match the factor levels exactly.",
           call. = FALSE)
    }
    matrix_value <- matrix_value[, levels, drop = FALSE]
  } else {
    colnames(matrix_value) <- levels
  }
  if (any(!is.finite(matrix_value))) {
    stop("Contrast coefficients must be finite.", call. = FALSE)
  }
  if (any(abs(rowSums(matrix_value)) > sqrt(tol))) {
    stop("Every location contrast must sum to zero across factor levels.",
         call. = FALSE)
  }
  if (any(rowSums(abs(matrix_value)) <= tol)) {
    stop("Every contrast must contain a non-zero coefficient.",
         call. = FALSE)
  }
  if (is.null(rownames(matrix_value))) {
    rownames(matrix_value) <- paste0("contrast_", seq_len(nrow(matrix_value)))
  }
  if (any(!nzchar(rownames(matrix_value))) || anyDuplicated(rownames(matrix_value))) {
    stop("Contrast names must be non-empty and unique.", call. = FALSE)
  }
  matrix_value
}

.contrast_projector <- function(coefficients, group, counts, H_full, tol) {
  score <- coefficients[as.character(group)] /
    as.numeric(counts[as.character(group)])
  denominator <- sum(score ^ 2)
  if (!is.finite(denominator) || denominator <= tol) {
    stop("A contrast has numerically zero estimable length.", call. = FALSE)
  }
  H <- tcrossprod(score) / denominator
  H_reduced <- H_full - H
  error <- max(
    abs(H %*% H - H),
    abs(H_reduced %*% H_reduced - H_reduced),
    abs(H %*% H_reduced)
  )
  if (error > max(100 * tol, sqrt(.Machine$double.eps))) {
    stop("A requested contrast does not form an estimable one-df projector.",
         call. = FALSE)
  }
  list(H = H, H_reduced = H_reduced, score = score, error = error)
}

#' Planned one-factor distance contrasts
#'
#' Tests one or more contrasts among the centroids of a single categorical
#' fixed factor. Each row supplies coefficients applied to unweighted group
#' centroids and must sum to zero. A contrast-specific reduced model retains all
#' orthogonal group information; Freedman--Lane residual permutations are
#' unrestricted. Holm adjustment is the default and remains valid for
#' non-orthogonal contrast families.
#'
#' @param design Canonical fixed one-factor Design Object.
#' @param contrasts Numeric contrast vector, matrix with contrasts in rows, or
#'   named list of numeric vectors.
#' @param adjust Multiplicity method: `"holm"`, `"bonferroni"`, `"BH"`,
#'   `"BY"`, or `"none"`.
#' @param permutations A non-negative integer or `"all"`.
#' @param seed Optional random seed.
#' @param keep_permutations Retain permutation statistics.
#' @param max_exact Maximum exact permutation space.
#' @param tol Numerical tolerance.
#' @return An object of class `permanova_contrasts`.
#' @export
permanova_contrasts <- function(
    design, contrasts, adjust = c("holm", "bonferroni", "BH", "BY", "none"),
    permutations = 999L, seed = NULL, keep_permutations = FALSE,
    max_exact = 100000L, tol = 1e-10) {
  adjust <- match.arg(adjust)
  structure <- .validate_contrast_design(design, tol = tol)
  contrast_matrix <- .as_contrast_matrix(
    contrasts, structure$levels, tol = tol
  )
  if (!is.null(seed)) set.seed(seed)
  plan <- .permutation_orders(design$n, permutations, max_exact)
  specifications <- lapply(seq_len(nrow(contrast_matrix)), function(index) {
    coefficients <- contrast_matrix[index, ]
    projector <- .contrast_projector(
      coefficients, structure$group, structure$counts,
      structure$H_full, tol
    )
    result <- .unit_freedman_lane(
      design$G, projector$H, 1L, projector$H_reduced,
      structure$H_denominator, structure$denominator_df, plan, tol
    )
    notes <- character()
    if (!isTRUE(design$spectral$euclidean)) {
      notes <- c(notes, "Empirical inference retained the complete non-Euclidean Gower spectrum.")
    }
    test <- list(
      term = structure$fixed$name,
      contrast = rownames(contrast_matrix)[[index]],
      coefficients = coefficients,
      hypothesis = paste0(
        "the centroid contrast `", rownames(contrast_matrix)[[index]],
        "` is zero"
      ),
      numerator = list(
        ss = result$numerator_ss, df = 1L, ms = result$numerator_ms,
        projector = projector$H
      ),
      denominator = list(
        expression = "MS[within-group Residual]",
        ss = result$denominator_ss, df = structure$denominator_df,
        ms = result$denominator_ms, projector = structure$H_denominator
      ),
      nuisance_model = "intercept plus every group-centroid direction orthogonal to the requested contrast",
      statistic = result$observed, p_value = result$p_value,
      p_adjusted = NA_real_, adjustment = adjust,
      exchangeable_units = "observations",
      permutation_strata = "unrestricted",
      permutation_algorithm = "contrast-specific Freedman-Lane residual permutation",
      attainable_permutations = plan$attainable,
      permutations_used = nrow(plan$orders),
      unique_permutations_used = if (nrow(plan$orders)) nrow(unique(plan$orders)) else 0L,
      exact_permutation = plan$exact && nrow(plan$orders) > 0L,
      approximation_status = if (!nrow(plan$orders)) "not tested" else if (plan$exact) {
        "exhaustive contrast-specific residual permutation"
      } else "Monte Carlo contrast-specific residual permutation",
      projector_error = projector$error, warnings = notes,
      permutation_statistics = if (keep_permutations) result$statistics else NULL
    )
    class(test) <- "permanova_test_spec"
    test
  })
  names(specifications) <- rownames(contrast_matrix)
  raw_p <- vapply(specifications, `[[`, numeric(1), "p_value")
  adjusted_p <- if (all(is.na(raw_p))) raw_p else
    stats::p.adjust(raw_p, method = adjust)
  for (index in seq_along(specifications)) {
    specifications[[index]]$p_adjusted <- adjusted_p[[index]]
  }
  table <- data.frame(
    contrast = names(specifications),
    df = 1L,
    ss = vapply(specifications, function(x) x$numerator$ss, numeric(1)),
    ms = vapply(specifications, function(x) x$numerator$ms, numeric(1)),
    denominator_df = structure$denominator_df,
    pseudo_f = vapply(specifications, `[[`, numeric(1), "statistic"),
    p_perm = raw_p, p_adjusted = adjusted_p,
    adjustment = adjust, stringsAsFactors = FALSE
  )
  out <- list(
    call = match.call(), design = design, factor = structure$fixed$name,
    levels = structure$levels, counts = structure$counts,
    contrast_matrix = contrast_matrix, table = table,
    tests = specifications, adjustment = adjust,
    family = "user-specified one-factor centroid contrasts",
    supported_design = paste(
      "one categorical fixed factor with positive within-group residual df;",
      "contrast-specific empirical inference"
    ),
    permutation = list(
      algorithm = "contrast-specific Freedman-Lane unrestricted residual permutation",
      exact = plan$exact && nrow(plan$orders) > 0L,
      attainable = plan$attainable, used = nrow(plan$orders), seed = seed
    )
  )
  class(out) <- "permanova_contrasts"
  out
}

#' Pairwise one-factor distance comparisons
#'
#' Constructs every pairwise difference among the levels of the sole fixed
#' factor and delegates to `permanova_contrasts()`.
#'
#' @inheritParams permanova_contrasts
#' @return An object of class `permanova_contrasts` with a pairwise family.
#' @export
permanova_pairwise <- function(
    design, adjust = c("holm", "bonferroni", "BH", "BY", "none"),
    permutations = 999L, seed = NULL, keep_permutations = FALSE,
    max_exact = 100000L, tol = 1e-10) {
  structure <- .validate_contrast_design(design, tol = tol)
  pairs <- utils::combn(structure$levels, 2L, simplify = FALSE)
  contrasts <- matrix(
    0, nrow = length(pairs), ncol = length(structure$levels),
    dimnames = list(NULL, structure$levels)
  )
  for (index in seq_along(pairs)) {
    contrasts[index, pairs[[index]][[1L]]] <- 1
    contrasts[index, pairs[[index]][[2L]]] <- -1
  }
  rownames(contrasts) <- vapply(
    pairs, function(pair) paste(pair, collapse = " - "), character(1)
  )
  out <- permanova_contrasts(
    design, contrasts = contrasts, adjust = match.arg(adjust),
    permutations = permutations, seed = seed,
    keep_permutations = keep_permutations,
    max_exact = max_exact, tol = tol
  )
  out$call <- match.call()
  out$family <- "all pairwise one-factor centroid differences"
  out
}

#' @export
print.permanova_contrasts <- function(x, ...) {
  cat("One-factor distance contrasts\n")
  cat("Family:", x$family, "\n")
  cat("Multiplicity adjustment:", x$adjustment, "\n")
  print(x$table, row.names = FALSE)
  invisible(x)
}
