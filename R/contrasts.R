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
