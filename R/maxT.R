.maxT_statistics <- function(design, structure, contrast_matrix, plan, tol) {
  projectors <- lapply(seq_len(nrow(contrast_matrix)), function(index) {
    .contrast_projector(
      contrast_matrix[index, ], structure$group, structure$counts,
      structure$H_full, tol
    )
  })
  denominator <- trace_hgh(structure$H_denominator, design$G) /
    structure$denominator_df
  if (!is.finite(denominator) || abs(denominator) <= tol) {
    stop("The within-group residual mean square is numerically zero.",
         call. = FALSE)
  }
  observed <- vapply(projectors, function(projector) {
    trace_hgh(projector$H, design$G) / denominator
  }, numeric(1))
  statistics <- matrix(
    numeric(nrow(plan$orders) * length(projectors)),
    nrow = nrow(plan$orders), ncol = length(projectors),
    dimnames = list(NULL, rownames(contrast_matrix))
  )
  if (nrow(plan$orders)) {
    progress <- .permutation_progress(
      nrow(plan$orders), "Permuting observations under the complete null"
    )
    on.exit(progress$close(), add = TRUE)
    for (row in seq_len(nrow(plan$orders))) {
      order <- plan$orders[row, ]
      G <- design$G[order, order, drop = FALSE]
      permuted_denominator <- trace_hgh(structure$H_denominator, G) /
        structure$denominator_df
      statistics[row, ] <- vapply(projectors, function(projector) {
        (trace_hgh(projector$H, G)) / permuted_denominator
      }, numeric(1))
      progress$tick(row)
    }
  }
  list(
    projectors = projectors, denominator = denominator,
    observed = observed, statistics = statistics,
    maxima = if (nrow(statistics)) apply(statistics, 1L, max) else numeric()
  )
}

#' Global-null max-F contrast inference
#'
#' Uses one shared raw-observation permutation distribution for every planned
#' contrast in a fixed one-factor design. Adjusted p-values compare each
#' observed pseudo-F with the permutation maximum. This controls familywise
#' error under the complete equal-centroid null (weak FWER); it does not claim
#' strong subset-null control.
#'
#' @param design Canonical fixed one-factor Design Object.
#' @param contrasts Numeric contrast vector, row matrix, or named list.
#' @param permutations A non-negative integer or `"all"`.
#' @param seed Optional random seed.
#' @param keep_permutations Retain the joint statistic matrix and maxima.
#' @param max_exact Maximum exact permutation space.
#' @param tol Numerical tolerance.
#' @return An object of class `permanova_maxT`.
#' @export
permanova_maxT <- function(
    design, contrasts, permutations = 999L, seed = NULL,
    keep_permutations = FALSE, max_exact = 100000L, tol = 1e-10) {
  structure <- .validate_contrast_design(design, tol = tol)
  contrast_matrix <- .as_contrast_matrix(contrasts, structure$levels, tol)
  if (!is.null(seed)) set.seed(seed)
  plan <- .permutation_orders(design$n, permutations, max_exact)
  joint <- .maxT_statistics(
    design, structure, contrast_matrix, plan, tol
  )
  raw <- adjusted <- rep(NA_real_, nrow(contrast_matrix))
  if (nrow(joint$statistics)) {
    for (index in seq_len(nrow(contrast_matrix))) {
      raw_extreme <- sum(
        joint$statistics[, index] >= joint$observed[[index]] - tol,
        na.rm = TRUE
      )
      max_extreme <- sum(
        joint$maxima >= joint$observed[[index]] - tol,
        na.rm = TRUE
      )
      if (plan$exact) {
        raw[[index]] <- raw_extreme / nrow(joint$statistics)
        adjusted[[index]] <- max_extreme / nrow(joint$statistics)
      } else {
        raw[[index]] <- (raw_extreme + 1) / (nrow(joint$statistics) + 1)
        adjusted[[index]] <- (max_extreme + 1) /
          (nrow(joint$statistics) + 1)
      }
    }
  }
  notes <- "max-F adjustment controls weak FWER under the complete equal-centroid null; strong subset-null control is not claimed."
  if (!isTRUE(design$spectral$euclidean)) {
    notes <- c(notes,
               "Empirical inference retained the complete non-Euclidean Gower spectrum.")
  }
  tests <- lapply(seq_len(nrow(contrast_matrix)), function(index) {
    projector <- joint$projectors[[index]]
    numerator_ss <- trace_hgh(projector$H, design$G)
    test <- list(
      term = structure$fixed$name,
      contrast = rownames(contrast_matrix)[[index]],
      coefficients = contrast_matrix[index, ],
      hypothesis = paste0("the centroid contrast `",
                          rownames(contrast_matrix)[[index]], "` is zero"),
      family_null = "all factor-level centroids are equal",
      numerator = list(
        ss = numerator_ss, df = 1L, ms = numerator_ss,
        projector = projector$H
      ),
      denominator = list(
        expression = "MS[within-group Residual]",
        ss = joint$denominator * structure$denominator_df,
        df = structure$denominator_df, ms = joint$denominator,
        projector = structure$H_denominator
      ),
      statistic = joint$observed[[index]], p_value = raw[[index]],
      p_adjusted = adjusted[[index]], adjustment = "global-null max-F",
      error_control = "weak familywise error control under the complete null",
      exchangeable_units = "observations",
      permutation_strata = "unrestricted",
      permutation_algorithm = "shared raw-observation permutation under the complete equal-centroid null",
      attainable_permutations = plan$attainable,
      permutations_used = nrow(plan$orders),
      exact_permutation = plan$exact && nrow(plan$orders) > 0L,
      approximation_status = if (!nrow(plan$orders)) {
        "not tested"
      } else if (plan$exact) {
        "exhaustive shared global-null permutation"
      } else "Monte Carlo shared global-null permutation",
      warnings = notes,
      permutation_statistics = if (keep_permutations) {
        joint$statistics[, index]
      } else NULL
    )
    class(test) <- "permanova_test_spec"
    test
  })
  names(tests) <- rownames(contrast_matrix)
  table <- data.frame(
    contrast = rownames(contrast_matrix), df = 1L,
    ss = vapply(tests, function(test) test$numerator$ss, numeric(1)),
    denominator_df = structure$denominator_df,
    pseudo_f = joint$observed, p_perm = raw, p_maxT = adjusted,
    error_control = "weak FWER under complete null",
    stringsAsFactors = FALSE
  )
  out <- list(
    call = match.call(), design = design, factor = structure$fixed$name,
    levels = structure$levels, counts = structure$counts,
    contrast_matrix = contrast_matrix, table = table, tests = tests,
    family = "user-specified one-factor centroid contrasts",
    family_null = "all factor-level centroids are equal",
    adjustment = "global-null max-F",
    error_control = "weak familywise error control under the complete equal-centroid null; strong control is not claimed",
    permutation = list(
      algorithm = "shared raw-observation global-null permutation",
      exact = plan$exact && nrow(plan$orders) > 0L,
      attainable = plan$attainable, used = nrow(plan$orders), seed = seed,
      statistic_matrix = if (keep_permutations) joint$statistics else NULL,
      maxima = if (keep_permutations) joint$maxima else NULL
    )
  )
  class(out) <- "permanova_maxT"
  out
}

#' Pairwise global-null max-F inference
#'
#' @inheritParams permanova_maxT
#' @return An object of class `permanova_maxT`.
#' @export
permanova_maxT_pairwise <- function(
    design, permutations = 999L, seed = NULL,
    keep_permutations = FALSE, max_exact = 100000L, tol = 1e-10) {
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
  out <- permanova_maxT(
    design, contrasts, permutations = permutations, seed = seed,
    keep_permutations = keep_permutations,
    max_exact = max_exact, tol = tol
  )
  out$call <- match.call()
  out$family <- "all pairwise one-factor centroid contrasts"
  out
}

#' @export
print.permanova_maxT <- function(x, ...) {
  cat("Global-null max-F distance contrasts\n")
  cat("Family:", x$family, "\n")
  cat("Error control:", x$error_control, "\n")
  print(x$table, row.names = FALSE)
  invisible(x)
}
