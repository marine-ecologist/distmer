.mixed_maxT_statistics <- function(structure, units, contrast_matrix, plan,
                                   tol) {
  projectors <- lapply(seq_len(nrow(contrast_matrix)), function(index) {
    coefficient_row <- as.numeric(
      contrast_matrix[index, ] %*% units$mean_rows
    )
    .linear_contrast_projector(
      units$X, coefficient_row, units$H_full, tol
    )
  })
  G <- structure$aggregated$G
  denominator <- trace_hgh(units$H_denominator, G) / units$denominator_df
  if (.mean_square_is_zero(denominator, G, tol)) {
    stop("The experimental-unit residual mean square is numerically zero.",
         call. = FALSE)
  }
  observed <- vapply(projectors, function(projector) {
    trace_hgh(projector$H, G) / denominator
  }, numeric(1))
  statistics <- matrix(
    numeric(nrow(plan$orders) * length(projectors)),
    nrow = nrow(plan$orders), ncol = length(projectors),
    dimnames = list(NULL, rownames(contrast_matrix))
  )
  if (nrow(plan$orders)) {
    progress <- .permutation_progress(
      nrow(plan$orders), "Permuting experimental units under the complete null"
    )
    on.exit(progress$close(), add = TRUE)
    for (row in seq_len(nrow(plan$orders))) {
      order <- plan$orders[row, ]
      permuted <- G[order, order, drop = FALSE]
      permuted_denominator <- trace_hgh(
        units$H_denominator, permuted
      ) / units$denominator_df
      if (.mean_square_is_zero(permuted_denominator, permuted, tol)) {
        stop("A permuted experimental-unit residual mean square is numerically zero.",
             call. = FALSE)
      }
      statistics[row, ] <- vapply(projectors, function(projector) {
        trace_hgh(projector$H, permuted) / permuted_denominator
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

#' Global-null max-F inference for supported mixed designs
#'
#' Uses one shared raw experimental-unit permutation distribution for a family
#' of treatment contrasts. Units and strata are selected from the declared
#' random structure. Adjusted p-values compare each observed pseudo-F with the
#' permutation maximum. This controls familywise error under the complete
#' treatment-null randomization hypothesis (weak FWER); strong subset-null
#' control is not claimed.
#'
#' @param design Canonical mixed Design Object containing a distance.
#' @param contrasts Numeric treatment-level contrast vector, row matrix, or
#'   named list.
#' @param permutations A non-negative integer or `"all"`.
#' @param seed Optional random seed.
#' @param keep_permutations Retain the joint statistic matrix and maxima.
#' @param max_exact Maximum exact permutation space.
#' @param tol Numerical tolerance.
#' @return An object of class `permanova_mixed_maxT`.
#' @export
permanova_mixed_maxT <- function(
    design, contrasts, permutations = 999L, seed = NULL,
    keep_permutations = FALSE, max_exact = 100000L, tol = 1e-10) {
  structure <- .mixed_contrast_structure(design, tol = tol)
  units <- .prepare_mixed_contrast_units(structure, tol = tol)
  contrast_matrix <- .as_contrast_matrix(contrasts, units$levels, tol)
  if (!is.null(seed)) set.seed(seed)
  plan <- .mixed_contrast_plan(structure, permutations, max_exact)
  joint <- .mixed_maxT_statistics(
    structure, units, contrast_matrix, plan, tol
  )
  raw <- adjusted <- rep(NA_real_, nrow(contrast_matrix))
  if (nrow(joint$statistics)) {
    for (index in seq_len(nrow(contrast_matrix))) {
      raw_extreme <- sum(
        joint$statistics[, index] >=
          joint$observed[[index]] - max(100 * tol, sqrt(.Machine$double.eps)) *
            max(1, abs(joint$observed[[index]]))
      )
      max_extreme <- sum(
        joint$maxima >=
          joint$observed[[index]] - max(100 * tol, sqrt(.Machine$double.eps)) *
            max(1, abs(joint$observed[[index]]))
      )
      if (plan$exact) {
        raw[[index]] <- raw_extreme / nrow(joint$statistics)
        adjusted[[index]] <- max_extreme / nrow(joint$statistics)
      } else {
        raw[[index]] <- (raw_extreme + 1) /
          (nrow(joint$statistics) + 1)
        adjusted[[index]] <- (max_extreme + 1) /
          (nrow(joint$statistics) + 1)
      }
    }
  }
  notes <- paste(
    "max-F adjustment controls weak FWER under the complete treatment-null",
    "randomization hypothesis; strong subset-null control is not claimed."
  )
  if (is.finite(plan$attainable) && plan$attainable < 100) {
    notes <- c(notes, paste0(
      "Only ", plan$attainable,
      " experimental-unit permutations are attainable."
    ))
  }
  if (!isTRUE(design$spectral$euclidean)) {
    notes <- c(
      notes,
      "Empirical inference retained the complete non-Euclidean Gower spectrum."
    )
  }
  strata <- if (is.null(structure$strata_name)) {
    "unrestricted"
  } else structure$strata_name
  algorithm <- paste0(
    "shared raw permutation of ", structure$unit_name,
    if (identical(strata, "unrestricted")) {
      " under the complete treatment null"
    } else paste0(" within ", strata, " under the complete treatment null")
  )
  tests <- lapply(seq_len(nrow(contrast_matrix)), function(index) {
    projector <- joint$projectors[[index]]
    numerator_ss <- trace_hgh(projector$H, structure$aggregated$G)
    test <- list(
      term = structure$fixed$name,
      contrast = rownames(contrast_matrix)[[index]],
      coefficients = contrast_matrix[index, ],
      hypothesis = paste0("the experimental-unit centroid contrast `",
                          rownames(contrast_matrix)[[index]], "` is zero"),
      family_null = "all treatment experimental-unit centroids are equal after declared model-stratum adjustment",
      estimand = structure$estimand,
      numerator = list(
        ss = numerator_ss, df = 1L, ms = numerator_ss,
        projector = projector$H
      ),
      denominator = list(
        expression = paste0("MS[", structure$denominator_name, "]"),
        ss = joint$denominator * units$denominator_df,
        df = units$denominator_df, ms = joint$denominator,
        projector = units$H_denominator
      ),
      statistic = joint$observed[[index]], p_value = raw[[index]],
      p_adjusted = adjusted[[index]], adjustment = "global-null max-F",
      error_control = "weak familywise error control under the complete treatment null",
      design_type = structure$type,
      exchangeable_units = structure$unit_name,
      permutation_strata = strata,
      permutation_algorithm = algorithm,
      attainable_permutations = plan$attainable,
      permutations_used = nrow(plan$orders),
      exact_permutation = plan$exact && nrow(plan$orders) > 0L,
      approximation_status = if (!nrow(plan$orders)) {
        "not tested"
      } else if (plan$exact) {
        "exhaustive shared experimental-unit global-null permutation"
      } else "Monte Carlo shared experimental-unit global-null permutation",
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
    contrast = rownames(contrast_matrix), estimand = structure$estimand,
    df = 1L,
    ss = vapply(tests, function(test) test$numerator$ss, numeric(1)),
    denominator = paste0("MS[", structure$denominator_name, "]"),
    denominator_df = units$denominator_df,
    pseudo_f = joint$observed, p_perm = raw, p_maxT = adjusted,
    error_control = "weak FWER under complete treatment null",
    stringsAsFactors = FALSE
  )
  out <- list(
    call = match.call(), design = design, design_type = structure$type,
    factor = structure$fixed$name, levels = units$levels,
    contrast_matrix = contrast_matrix, table = table, tests = tests,
    family = "user-specified mixed-design experimental-unit contrasts",
    family_null = "all treatment experimental-unit centroids are equal after declared model-stratum adjustment",
    estimand = structure$estimand,
    denominator = paste0("MS[", structure$denominator_name, "]"),
    adjustment = "global-null max-F",
    error_control = paste(
      "weak familywise error control under the complete treatment-null",
      "randomization hypothesis; strong control is not claimed"
    ),
    permutation = list(
      algorithm = algorithm, units = structure$unit_name,
      strata = structure$strata_name,
      exact = plan$exact && nrow(plan$orders) > 0L,
      attainable = plan$attainable, used = nrow(plan$orders), seed = seed,
      statistic_matrix = if (keep_permutations) joint$statistics else NULL,
      maxima = if (keep_permutations) joint$maxima else NULL
    )
  )
  class(out) <- "permanova_mixed_maxT"
  out
}

#' Pairwise global-null max-F inference for supported mixed designs
#'
#' @inheritParams permanova_mixed_maxT
#' @return An object of class `permanova_mixed_maxT`.
#' @export
permanova_mixed_maxT_pairwise <- function(
    design, permutations = 999L, seed = NULL,
    keep_permutations = FALSE, max_exact = 100000L, tol = 1e-10) {
  structure <- .mixed_contrast_structure(design, tol = tol)
  levels <- levels(droplevels(structure$treatment))
  pairs <- utils::combn(levels, 2L, simplify = FALSE)
  contrasts <- matrix(
    0, nrow = length(pairs), ncol = length(levels),
    dimnames = list(NULL, levels)
  )
  for (index in seq_along(pairs)) {
    contrasts[index, pairs[[index]][[1L]]] <- 1
    contrasts[index, pairs[[index]][[2L]]] <- -1
  }
  rownames(contrasts) <- vapply(
    pairs, function(pair) paste(pair, collapse = " - "), character(1)
  )
  out <- permanova_mixed_maxT(
    design, contrasts, permutations = permutations, seed = seed,
    keep_permutations = keep_permutations,
    max_exact = max_exact, tol = tol
  )
  out$call <- match.call()
  out$family <- "all pairwise mixed-design experimental-unit contrasts"
  out
}

#' @export
print.permanova_mixed_maxT <- function(x, ...) {
  cat("Mixed-design global-null max-F distance contrasts\n")
  cat("Design:", x$design_type, "\n")
  cat("Estimand:", x$estimand, "\n")
  cat("Error control:", x$error_control, "\n")
  print(x$table, row.names = FALSE)
  invisible(x)
}
