.unit_factor_from_blocks <- function(value, blocks, label) {
  labels <- vapply(blocks, function(index) {
    observed <- unique(as.character(value[index]))
    if (length(observed) != 1L) {
      stop("`", label, "` must be constant within every contrast unit.",
           call. = FALSE)
    }
    observed
  }, character(1))
  factor(labels, levels = levels(droplevels(value)))
}

.mixed_contrast_structure <- function(design, tol = 1e-10) {
  if (!inherits(design, "permanova_design") || is.null(design$G)) {
    stop("`design` must be a compiled Design Object containing `distance`.",
         call. = FALSE)
  }
  fixed <- .one_categorical_fixed(design, "Mixed-design contrasts")
  origins <- if (length(design$random_terms)) {
    vapply(design$random_terms, `[[`, character(1), "origin")
  } else character()

  if (length(origins) == 1L && identical(origins, "main")) {
    cluster_label <- design$random_terms[[1L]]$label
    cluster_sizes <- as.integer(table(design$groups[[cluster_label]]))
    if (length(unique(cluster_sizes)) != 1L) {
      validated <- .validate_unequal_random_intercept(design, tol = tol)
      return(list(
        type = "unequally replicated random intercept", fixed = fixed,
        validated = validated, aggregated = validated$aggregated,
        treatment = validated$treatment,
        model_strata = NULL,
        permutation_strata = validated$size_strata,
        strata_name = "exact random-intercept unit observation count",
        unit_name = validated$cluster_name,
        estimand = "equally weighted random-intercept unit centroids",
        denominator_name = "unit-centroid Residual",
        permutation_algorithm = "contrast-specific Freedman-Lane permutation of whole random-intercept units within exact size strata"
      ))
    }
    validated <- .validate_balanced_random_intercept(design, tol = tol)
    aggregated <- .aggregate_gower_units(design, validated$cluster, tol = tol)
    treatment <- .unit_factor_from_blocks(
      fixed$value, aggregated$blocks, fixed$name
    )
    return(list(
      type = "balanced random intercept", fixed = fixed,
      validated = validated, aggregated = aggregated,
      treatment = treatment, model_strata = NULL,
      permutation_strata = NULL, strata_name = NULL,
      unit_name = validated$cluster_name,
      estimand = "equally weighted random-intercept unit centroids",
      denominator_name = paste0(validated$cluster_name,
                                " within treatment model"),
      permutation_algorithm = "contrast-specific Freedman-Lane permutation of whole random-intercept units"
    ))
  }

  if (length(origins) >= 2L && origins[[1L]] == "main" &&
      all(origins[-1L] == "nesting")) {
    hierarchy_labels <- vapply(
      design$random_terms, `[[`, character(1), "label"
    )
    if (!all(design$balance$random_balanced[hierarchy_labels])) {
      if (length(origins) != 2L) {
        stop("Unbalanced mixed contrasts support exactly one parent/child hierarchy.",
             call. = FALSE)
      }
      assignment <- .unbalanced_nested_assignment(design)$level
      if (assignment != "top") {
        lower <- .validate_unbalanced_nested_lower(
          design, level = assignment, tol = tol
        )
        return(lower$structure)
      }
      validated <- .validate_unbalanced_nested(design, tol = tol)
      return(list(
        type = "unbalanced top-level parent/child hierarchy",
        fixed = fixed, validated = validated,
        aggregated = validated$aggregated,
        treatment = validated$treatment,
        model_strata = NULL,
        permutation_strata = validated$topology_strata,
        strata_name = "identical sorted child observation-count topology",
        unit_name = validated$parent_name,
        estimand = "equally weighted intact top-level hierarchy-unit centroids",
        denominator_name = "top-unit-centroid Residual",
        permutation_algorithm = "contrast-specific Freedman-Lane permutation of intact top-level units within identical hierarchy topologies"
      ))
    }
    validated <- .validate_nested_design(design, tol = tol)
    if (validated$term_levels[[fixed$name]] != 1L) {
      stop("Nested mixed contrasts require treatment to lie at the top hierarchy level.",
           call. = FALSE)
    }
    aggregated <- .aggregate_gower_units(design, validated$top, tol = tol)
    treatment <- .unit_factor_from_blocks(
      fixed$value, aggregated$blocks, fixed$name
    )
    return(list(
      type = "balanced top-level nested hierarchy", fixed = fixed,
      validated = validated, aggregated = aggregated,
      treatment = treatment, model_strata = NULL,
      permutation_strata = NULL, strata_name = NULL,
      unit_name = validated$top_name,
      estimand = "equally weighted intact top-level hierarchy-unit centroids",
      denominator_name = validated$denominator_name,
      permutation_algorithm = "contrast-specific Freedman-Lane permutation of intact top-level hierarchy units"
    ))
  }

  if (length(origins) == 2L && sum(origins == "main") == 1L &&
      sum(origins == "interaction") == 1L) {
    main_term <- design$random_terms[[which(origins == "main")]]
    interaction_term <- design$random_terms[[which(origins == "interaction")]]
    block_factor <- droplevels(design$groups[[main_term$label]])
    cell_factor <- droplevels(design$groups[[interaction_term$label]])
    mapping <- unique(data.frame(
      cell = as.integer(cell_factor), block = as.integer(block_factor),
      treatment = as.integer(fixed$value)
    ))
    incidence <- table(mapping$block, mapping$treatment)
    complete <- nrow(mapping) == nlevels(block_factor) * nlevels(fixed$value) &&
      length(incidence) == nlevels(block_factor) * nlevels(fixed$value) &&
      all(incidence == 1L)
    equal_cells <- length(unique(as.integer(table(cell_factor)))) == 1L
    if (!complete || !equal_cells) {
      validated <- .validate_unbalanced_block(design, tol = tol)
      return(list(
        type = "incomplete or unequally replicated randomized block",
        fixed = fixed, validated = validated,
        aggregated = validated$aggregated,
        treatment = validated$treatment,
        model_strata = validated$block_factor,
        permutation_strata = validated$exchangeability_strata,
        strata_name = "block crossed with exact treatment-cell size",
        unit_name = validated$cell_name,
        estimand = "equally weighted observed treatment-by-block cell centroids with additive block adjustment",
        denominator_name = "cell-centroid Residual after block and treatment",
        permutation_algorithm = "contrast-specific Freedman-Lane permutation of cell-centroid residuals within block-by-size strata"
      ))
    }
    crossed <- .validate_crossed_design(
      design, assume_additive = FALSE, tol = tol
    )
    validated <- .validate_crossed_permutation_design(
      design, crossed, assume_additive = FALSE, tol = tol
    )
    aggregated <- .aggregate_gower_units(design, validated$cell, tol = tol)
    treatment <- .unit_factor_from_blocks(
      fixed$value, aggregated$blocks, fixed$name
    )
    block <- .unit_factor_from_blocks(
      design$data[[validated$block_name]], aggregated$blocks,
      validated$block_name
    )
    return(list(
      type = "balanced treatment by random block", fixed = fixed,
      validated = validated, aggregated = aggregated,
      treatment = treatment, model_strata = block,
      permutation_strata = block,
      strata_name = validated$block_name,
      unit_name = validated$interaction_name,
      estimand = "equally weighted treatment-by-block cell centroids",
      denominator_name = paste0(fixed$name, ":", validated$block_name),
      permutation_algorithm = "contrast-specific Freedman-Lane permutation of treatment-cell centroids within random blocks"
    ))
  }

  if (length(origins) == 4L && sum(origins == "main") == 1L &&
      sum(origins == "nesting") == 1L &&
      sum(origins == "interaction") == 2L) {
    validated <- tryCatch(
      .validate_hybrid_design(design, tol = tol),
      error = function(error) error
    )
    if (inherits(validated, "error")) {
      return(.validate_unbalanced_hybrid(design, tol = tol)$structure)
    }
    aggregated <- .aggregate_gower_units(
      design, validated$treatment_region, tol = tol
    )
    treatment <- .unit_factor_from_blocks(
      fixed$value, aggregated$blocks, fixed$name
    )
    region <- .unit_factor_from_blocks(
      design$data[[validated$region_name]], aggregated$blocks,
      validated$region_name
    )
    return(list(
      type = "balanced crossed-nested treatment by region/site hybrid",
      fixed = fixed, validated = validated, aggregated = aggregated,
      treatment = treatment, model_strata = region,
      permutation_strata = region,
      strata_name = validated$region_name,
      unit_name = validated$treatment_region_term,
      estimand = "equally weighted treatment-by-region profile centroids",
      denominator_name = "treatment:region",
      permutation_algorithm = "contrast-specific Freedman-Lane permutation of treatment profiles within regions"
    ))
  }

  stop(
    "Mixed-design contrasts support one random intercept, one parent/child hierarchy, a declared treatment-by-block design, or a validated crossed-nested hybrid.",
    call. = FALSE
  )
}

.prepare_mixed_contrast_units <- function(structure, tol = 1e-10) {
  unit_data <- data.frame(.treatment = structure$treatment)
  if (!is.null(structure$model_strata)) {
    unit_data$.stratum <- structure$model_strata
  }
  formula <- if (is.null(structure$model_strata)) {
    ~ .treatment
  } else {
    ~ .stratum + .treatment
  }
  X <- stats::model.matrix(formula, unit_data)
  full <- projection_matrix(X, tol = tol)
  if (full$rank != ncol(X)) {
    stop("The experimental-unit treatment model is rank deficient.",
         call. = FALSE)
  }
  denominator_df <- nrow(X) - full$rank
  if (denominator_df <= 0L) {
    stop("No experimental-unit residual degrees of freedom remain for mixed contrasts.",
         call. = FALSE)
  }
  levels <- levels(droplevels(structure$treatment))
  if (is.null(structure$model_strata)) {
    mean_rows <- do.call(rbind, lapply(levels, function(level) {
      colMeans(X[structure$treatment == level, , drop = FALSE])
    }))
  } else {
    reference <- expand.grid(
      .treatment = levels,
      .stratum = levels(droplevels(structure$model_strata)),
      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
    )
    reference$.treatment <- factor(reference$.treatment, levels = levels)
    reference$.stratum <- factor(
      reference$.stratum,
      levels = levels(droplevels(structure$model_strata))
    )
    reference_X <- stats::model.matrix(formula, reference)
    if (!identical(colnames(reference_X), colnames(X))) {
      stop("Experimental-unit reference columns do not match the fitted model.",
           call. = FALSE)
    }
    mean_rows <- do.call(rbind, lapply(levels, function(level) {
      colMeans(reference_X[reference$.treatment == level, , drop = FALSE])
    }))
  }
  rownames(mean_rows) <- levels
  colnames(mean_rows) <- colnames(X)
  list(
    data = unit_data, formula = formula, X = X, H_full = full$H,
    H_denominator = diag(nrow(X)) - full$H,
    denominator_df = denominator_df, levels = levels,
    mean_rows = mean_rows
  )
}

.mixed_contrast_plan <- function(structure, permutations, max_exact) {
  if (is.null(structure$permutation_strata)) {
    .permutation_orders(
      length(structure$treatment), permutations = permutations,
      max_exact = max_exact
    )
  } else {
    .stratified_permutation_plan(
      structure$permutation_strata, permutations = permutations,
      max_exact = max_exact
    )
  }
}

#' Contrasts for supported mixed PERMANOVA designs
#'
#' Tests planned contrasts among a single categorical treatment's experimental-
#' unit centroids. The unit and exact interaction/residual denominator are
#' selected from the declared random structure. Supported structures are one
#' balanced or unequally replicated random intercept, a top-level treatment in
#' one balanced or unbalanced parent/child hierarchy, a complete or incomplete
#' treatment-by-block design with its declared random interaction, and the
#' balanced or validated unbalanced treatment-by-region/site hybrid.
#'
#' @param design Canonical mixed Design Object containing a distance.
#' @param contrasts Numeric treatment-level contrast vector, row matrix, or
#'   named list.
#' @param adjust Multiplicity adjustment.
#' @param permutations A non-negative integer or `"all"`.
#' @param seed Optional random seed.
#' @param keep_permutations Retain permutation statistics.
#' @param max_exact Maximum exact permutation space.
#' @param tol Numerical tolerance.
#' @return An object of class `permanova_mixed_contrasts`.
#' @export
permanova_mixed_contrasts <- function(
    design, contrasts,
    adjust = c("holm", "bonferroni", "BH", "BY", "none"),
    permutations = 999L, seed = NULL, keep_permutations = FALSE,
    max_exact = 100000L, tol = 1e-10) {
  adjust <- match.arg(adjust)
  structure <- .mixed_contrast_structure(design, tol = tol)
  units <- .prepare_mixed_contrast_units(structure, tol = tol)
  contrast_matrix <- .as_contrast_matrix(contrasts, units$levels, tol)
  if (!is.null(seed)) set.seed(seed)
  plan <- .mixed_contrast_plan(structure, permutations, max_exact)

  tests <- lapply(seq_len(nrow(contrast_matrix)), function(index) {
    coefficients <- contrast_matrix[index, ]
    coefficient_row <- as.numeric(coefficients %*% units$mean_rows)
    projector <- .linear_contrast_projector(
      units$X, coefficient_row, units$H_full, tol
    )
    result <- .unit_freedman_lane(
      structure$aggregated$G, projector$H, 1L, projector$H_reduced,
      units$H_denominator, units$denominator_df, plan, tol,
      label = paste0("Permuting unit contrast `",
                     rownames(contrast_matrix)[[index]], "`")
    )
    notes <- result$resolution_note
    if (is.finite(plan$attainable) && plan$attainable < 100) {
      notes <- c(notes, paste0(
        "Only ", plan$attainable,
        " experimental-unit permutations are attainable."
      ))
    }
    if (!isTRUE(design$spectral$euclidean)) {
      notes <- c(notes,
                 "Empirical inference retained the complete non-Euclidean Gower spectrum.")
    }
    test <- list(
      term = structure$fixed$name,
      contrast = rownames(contrast_matrix)[[index]],
      coefficients = coefficients,
      hypothesis = paste0("the experimental-unit centroid contrast `",
                          rownames(contrast_matrix)[[index]], "` is zero"),
      estimand = structure$estimand,
      numerator = list(
        ss = result$numerator_ss, df = 1L, ms = result$numerator_ms,
        projector = projector$H
      ),
      denominator = list(
        expression = paste0("MS[", structure$denominator_name, "]"),
        ss = result$denominator_ss, df = units$denominator_df,
        ms = result$denominator_ms, projector = units$H_denominator
      ),
      statistic = result$observed, p_value = result$p_value,
      p_adjusted = NA_real_, adjustment = adjust,
      design_type = structure$type,
      exchangeable_units = structure$unit_name,
      permutation_strata = if (is.null(structure$strata_name)) {
        "unrestricted"
      } else structure$strata_name,
      permutation_algorithm = structure$permutation_algorithm,
      attainable_permutations = plan$attainable,
      permutations_used = nrow(plan$orders),
      exact_permutation = plan$exact && nrow(plan$orders) > 0L,
      approximation_status = if (!nrow(plan$orders)) {
        "not tested"
      } else if (plan$exact) {
        "exhaustive experimental-unit residual permutation"
      } else "Monte Carlo experimental-unit residual permutation",
      projector_error = projector$error, warnings = notes,
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
    contrast = names(tests), estimand = structure$estimand, df = 1L,
    ss = vapply(tests, function(x) x$numerator$ss, numeric(1)),
    denominator = paste0("MS[", structure$denominator_name, "]"),
    denominator_df = units$denominator_df,
    pseudo_f = vapply(tests, `[[`, numeric(1), "statistic"),
    p_perm = raw, p_adjusted = adjusted, adjustment = adjust,
    stringsAsFactors = FALSE
  )
  out <- list(
    call = match.call(), design = design, design_type = structure$type,
    factor = structure$fixed$name, levels = units$levels,
    contrast_matrix = contrast_matrix, table = table, tests = tests,
    family = "user-specified mixed-design experimental-unit contrasts",
    estimand = structure$estimand,
    denominator = paste0("MS[", structure$denominator_name, "]"),
    supported_design = structure$type,
    permutation = list(
      algorithm = structure$permutation_algorithm,
      units = structure$unit_name, strata = structure$strata_name,
      exact = plan$exact && nrow(plan$orders) > 0L,
      attainable = plan$attainable, used = nrow(plan$orders), seed = seed
    )
  )
  class(out) <- "permanova_mixed_contrasts"
  out
}

#' Pairwise contrasts for supported mixed PERMANOVA designs
#'
#' @inheritParams permanova_mixed_contrasts
#' @return An object of class `permanova_mixed_contrasts`.
#' @export
permanova_mixed_pairwise <- function(
    design, adjust = c("holm", "bonferroni", "BH", "BY", "none"),
    permutations = 999L, seed = NULL, keep_permutations = FALSE,
    max_exact = 100000L, tol = 1e-10) {
  structure <- .mixed_contrast_structure(design, tol = tol)
  levels <- levels(droplevels(structure$treatment))
  pairs <- utils::combn(levels, 2L, simplify = FALSE)
  contrast_matrix <- matrix(
    0, nrow = length(pairs), ncol = length(levels),
    dimnames = list(NULL, levels)
  )
  for (index in seq_along(pairs)) {
    contrast_matrix[index, pairs[[index]][[1L]]] <- 1
    contrast_matrix[index, pairs[[index]][[2L]]] <- -1
  }
  rownames(contrast_matrix) <- vapply(
    pairs, function(pair) paste(pair, collapse = " - "), character(1)
  )
  out <- permanova_mixed_contrasts(
    design, contrasts = contrast_matrix, adjust = match.arg(adjust),
    permutations = permutations, seed = seed,
    keep_permutations = keep_permutations,
    max_exact = max_exact, tol = tol
  )
  out$call <- match.call()
  out$family <- "all pairwise mixed-design experimental-unit contrasts"
  out
}

#' @export
print.permanova_mixed_contrasts <- function(x, ...) {
  cat("Mixed-design distance contrasts\n")
  cat("Design:", x$design_type, "\n")
  cat("Estimand:", x$estimand, "\n")
  cat("Denominator:", x$denominator, "\n")
  print(x$table, row.names = FALSE)
  invisible(x)
}
