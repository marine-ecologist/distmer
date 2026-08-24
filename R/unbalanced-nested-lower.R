.unbalanced_nested_assignment <- function(design) {
  if (!inherits(design, "permanova_design") || is.null(design$G)) {
    stop("`design` must be a compiled Design Object containing `distance`.",
         call. = FALSE)
  }
  fixed <- .one_categorical_fixed(design, "Lower-level unbalanced nested inference")
  if (length(design$random_terms) != 2L) {
    stop("Lower-level unbalanced nested inference requires exactly one parent/child hierarchy.",
         call. = FALSE)
  }
  origins <- vapply(design$random_terms, `[[`, character(1), "origin")
  parent_terms <- design$random_terms[origins == "main"]
  child_terms <- design$random_terms[origins == "nesting"]
  if (length(parent_terms) != 1L || length(child_terms) != 1L ||
      !identical(child_terms[[1L]]$parent, parent_terms[[1L]]$label)) {
    stop("Random terms must be one `/`-generated parent/child hierarchy.",
         call. = FALSE)
  }
  parent_name <- parent_terms[[1L]]$label
  child_name <- child_terms[[1L]]$label
  parent <- droplevels(design$groups[[parent_name]])
  child <- droplevels(design$groups[[child_name]])
  mapping <- unique(data.frame(
    child = as.integer(child), parent = as.integer(parent)
  ))
  if (any(table(mapping$child) != 1L)) {
    stop("Every nested child must belong to exactly one parent.", call. = FALSE)
  }
  constant_within <- function(group) {
    all(vapply(split(fixed$value, group), function(value) {
      length(unique(as.character(value))) == 1L
    }, logical(1)))
  }
  level <- if (constant_within(parent)) {
    "top"
  } else if (constant_within(child)) {
    "child"
  } else "observation"
  list(
    fixed = fixed, parent_name = parent_name, child_name = child_name,
    parent = parent, child = child, mapping = mapping, level = level
  )
}

.validate_unbalanced_nested_lower <- function(
    design, level = c("auto", "child", "observation"), tol = 1e-10) {
  level <- match.arg(level)
  hierarchy <- .unbalanced_nested_assignment(design)
  hierarchy_labels <- vapply(
    design$random_terms, `[[`, character(1), "label"
  )
  if (all(design$balance$random_balanced[hierarchy_labels])) {
    stop(
      "Lower-level empirical inference requires an unbalanced hierarchy; use `permanova_nested()` for the balanced nested EMS model.",
      call. = FALSE
    )
  }
  detected <- hierarchy$level
  if (detected == "top") {
    stop(
      "Treatment is assigned at the top hierarchy level; use `permanova_unbalanced_nested()` or top-level mixed contrasts.",
      call. = FALSE
    )
  }
  if (level != "auto" && level != detected) {
    stop("Requested lower hierarchy level `", level,
         "` does not match detected `", detected, "` assignment.",
         call. = FALSE)
  }
  level <- detected
  if (level == "child") {
    aggregated <- .aggregate_gower_units(design, hierarchy$child, tol = tol)
    treatment <- .unit_factor_from_blocks(
      hierarchy$fixed$value, aggregated$blocks, hierarchy$fixed$name
    )
    model_strata <- .unit_factor_from_blocks(
      hierarchy$parent, aggregated$blocks, hierarchy$parent_name
    )
    permutation_strata <- interaction(
      model_strata, factor(aggregated$sizes),
      drop = TRUE, lex.order = TRUE
    )
    strata_name <- paste0(
      hierarchy$parent_name, " crossed with exact ",
      hierarchy$child_name, " observation count"
    )
    unit_name <- hierarchy$child_name
    estimand <- paste0(
      "equally weighted `", hierarchy$child_name,
      "` centroids adjusted for `", hierarchy$parent_name, "`"
    )
    denominator_name <- paste0(
      hierarchy$child_name,
      "-centroid Residual after parent and treatment"
    )
    algorithm <- paste0(
      "Freedman-Lane permutation of intact `", hierarchy$child_name,
      "`-centroid residuals within parent-by-exact-size strata"
    )
  } else {
    observation <- factor(seq_len(design$n))
    aggregated <- .aggregate_gower_units(design, observation, tol = tol)
    treatment <- droplevels(hierarchy$fixed$value)
    model_strata <- hierarchy$child
    permutation_strata <- hierarchy$child
    strata_name <- hierarchy$child_name
    unit_name <- "observations"
    estimand <- paste0(
      "observation-level treatment location adjusted for `",
      hierarchy$child_name, "`"
    )
    denominator_name <- paste0(
      "Observation Residual after ", hierarchy$child_name,
      " and treatment"
    )
    algorithm <- paste0(
      "Freedman-Lane permutation of observation residuals within `",
      hierarchy$child_name, "`"
    )
  }
  structure <- list(
    type = paste0("unbalanced nested ", level, "-level treatment"),
    fixed = hierarchy$fixed, aggregated = aggregated,
    treatment = treatment, model_strata = model_strata,
    permutation_strata = permutation_strata,
    strata_name = strata_name, unit_name = unit_name,
    estimand = estimand, denominator_name = denominator_name,
    permutation_algorithm = algorithm
  )
  units <- tryCatch(
    .prepare_mixed_contrast_units(structure, tol = tol),
    error = function(error) {
      if (grepl("rank deficient", conditionMessage(error), fixed = TRUE)) {
        stop(
          "The lower-level treatment incidence is disconnected or aliases at least one treatment contrast.",
          call. = FALSE
        )
      }
      stop(conditionMessage(error), call. = FALSE)
    }
  )
  reduced_X <- stats::model.matrix(~ .stratum, units$data)
  reduced <- projection_matrix(reduced_X, tol = tol)
  numerator_df <- projection_matrix(units$X, tol = tol)$rank - reduced$rank
  expected_df <- nlevels(droplevels(treatment)) - 1L
  if (numerator_df != expected_df) {
    stop(
      "The lower-level treatment incidence is disconnected or aliases at least one treatment contrast.",
      call. = FALSE
    )
  }
  treatment_columns <- stats::model.matrix(~ 0 + treatment)
  exchangeable_rank <- .stratified_treatment_rank(
    treatment_columns, permutation_strata, tol
  )
  if (exchangeable_rank < expected_df) {
    stop(
      if (level == "child") {
        "Parent-by-child-size strata do not span every treatment contrast; intact child-centroid permutation cannot test the full hypothesis."
      } else {
        "Within-child observation permutations do not span every treatment contrast; observation-level inference is confounded with child."
      }, call. = FALSE
    )
  }
  structure$validated <- list(
    assignment_level = level, parent_name = hierarchy$parent_name,
    child_name = hierarchy$child_name,
    parent = hierarchy$parent, child = hierarchy$child,
    unit_sizes = aggregated$sizes,
    exchangeable_rank = exchangeable_rank
  )
  list(
    hierarchy = hierarchy, assignment_level = level,
    structure = structure, units = units,
    H_reduced = reduced$H,
    H_numerator = units$H_full - reduced$H,
    numerator_df = numerator_df
  )
}

#' Lower-level unbalanced nested PERMANOVA
#'
#' Tests one categorical treatment assigned either to child units or to
#' observations in one unbalanced parent/child hierarchy. Child-level tests use
#' equally weighted child centroids, adjust for parent, and permute intact
#' children only within parent-by-exact-size strata. Observation-level tests
#' adjust for child and permute observation residuals within child.
#'
#' @param design Canonical unbalanced parent/child Design Object containing a
#'   distance.
#' @param level `"auto"`, `"child"`, or `"observation"`.
#' @param permutations A non-negative integer or `"all"`.
#' @param seed Optional random seed.
#' @param keep_permutations Retain permutation statistics.
#' @param max_exact Maximum exact permutation space.
#' @param tol Numerical tolerance.
#' @return An object of class `permanova_unbalanced_nested_lower`.
#' @export
permanova_unbalanced_nested_lower <- function(
    design, level = c("auto", "child", "observation"),
    permutations = 999L, seed = NULL, keep_permutations = FALSE,
    max_exact = 100000L, tol = 1e-10) {
  validated <- .validate_unbalanced_nested_lower(
    design, level = match.arg(level), tol = tol
  )
  structure <- validated$structure
  units <- validated$units
  if (!is.null(seed)) set.seed(seed)
  plan <- .mixed_contrast_plan(structure, permutations, max_exact)
  result <- .unit_freedman_lane(
    structure$aggregated$G, validated$H_numerator,
    validated$numerator_df, validated$H_reduced,
    units$H_denominator, units$denominator_df, plan, tol,
    label = "Permuting lower-level unit centroids"
  )
  notes <- result$resolution_note
  if (is.finite(plan$attainable) && plan$attainable < 100) {
    notes <- c(notes, paste0(
      "Only ", plan$attainable,
      " lower-level permutations are attainable."
    ))
  }
  if (!isTRUE(design$spectral$euclidean)) {
    notes <- c(
      notes,
      "Empirical inference retained the complete non-Euclidean Gower spectrum."
    )
  }
  term <- structure$fixed$name
  test <- list(
    term = term,
    hypothesis = paste0(
      "the `", term, "` effect is zero at the unbalanced nested `",
      validated$assignment_level, "` level"
    ),
    assignment_level = validated$assignment_level,
    estimand = structure$estimand,
    numerator = list(
      ss = result$numerator_ss, df = validated$numerator_df,
      ms = result$numerator_ms, projector = validated$H_numerator
    ),
    denominator = list(
      expression = paste0("MS[", structure$denominator_name, "]"),
      ss = result$denominator_ss, df = units$denominator_df,
      ms = result$denominator_ms, projector = units$H_denominator
    ),
    statistic = result$observed, p_value = result$p_value,
    exchangeable_units = structure$unit_name,
    permutation_strata = structure$strata_name,
    permutation_algorithm = structure$permutation_algorithm,
    attainable_permutations = plan$attainable,
    permutations_used = nrow(plan$orders),
    exact_permutation = plan$exact && nrow(plan$orders) > 0L,
    approximation_status = if (!nrow(plan$orders)) {
      "not tested"
    } else if (plan$exact) {
      "exhaustive lower-level restricted residual permutation"
    } else "Monte Carlo lower-level restricted residual permutation",
    warnings = notes,
    permutation_statistics = if (keep_permutations) {
      result$statistics
    } else NULL
  )
  class(test) <- "permanova_test_spec"
  table <- data.frame(
    term = c(term, structure$denominator_name),
    assignment_level = c(validated$assignment_level, NA_character_),
    df = c(validated$numerator_df, units$denominator_df),
    ss = c(result$numerator_ss, result$denominator_ss),
    ms = c(result$numerator_ms, result$denominator_ms),
    pseudo_f = c(result$observed, NA_real_),
    p_perm = c(result$p_value, NA_real_),
    stringsAsFactors = FALSE
  )
  out <- list(
    call = match.call(), design = design, table = table, test = test,
    assignment_level = validated$assignment_level,
    structure = list(
      parent = validated$hierarchy$parent_name,
      child = validated$hierarchy$child_name,
      units = nrow(units$X), unit_sizes = structure$aggregated$sizes
    ),
    estimand = structure$estimand,
    supported_design = paste0(
      "one categorical ", validated$assignment_level,
      "-level treatment in one unbalanced parent/child hierarchy"
    ),
    permutation = list(
      algorithm = structure$permutation_algorithm,
      units = structure$unit_name, strata = structure$strata_name,
      exact = test$exact_permutation, attainable = plan$attainable,
      used = nrow(plan$orders), seed = seed
    )
  )
  class(out) <- "permanova_unbalanced_nested_lower"
  out
}

#' @export
print.permanova_unbalanced_nested_lower <- function(x, ...) {
  cat("Lower-level unbalanced nested PERMANOVA\n")
  cat("Assignment level:", x$assignment_level, "\n")
  cat("Estimand:", x$estimand, "\n")
  cat("Permutation:", x$permutation$algorithm, "\n")
  print(x$table, row.names = FALSE)
  invisible(x)
}
