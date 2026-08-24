# Multi-level unbalanced nested inference.
#
# The two-level engines (`permanova_unbalanced_nested`, and its lower-level
# companion) aggregate observations straight onto parent or child centroids. A
# chain deeper than parent/child cannot be collapsed that way: a single
# observation-weighted mean at the top of a `reg/site/quad` chain silently
# reweights regions by how many observations they happen to carry, which is the
# imbalance the equal-weight estimand exists to remove.
#
# Aggregation is therefore cascaded. Observations are averaged onto the deepest
# units, those unit centroids are averaged onto their parents with equal weight
# per child, and so on up to the level being tested. Because every stage is a
# linear map, the whole cascade composes into one membership matrix `M`, and the
# aggregated Gower matrix is `J M G M' J` exactly as in the two-level case.

# Ordered `/`-generated chain, top level first. Returns NULL when the design is
# not one contiguous nesting chain, so callers can fall through to their own
# error message rather than inheriting this one.
.nested_chain <- function(design) {
  if (!inherits(design, "permanova_design") || is.null(design$G)) return(NULL)
  terms <- design$random_terms
  if (length(terms) < 2L) return(NULL)
  origins <- vapply(terms, `[[`, character(1), "origin")
  if (origins[[1L]] != "main" || any(origins[-1L] != "nesting")) return(NULL)
  labels <- vapply(terms, `[[`, character(1), "label")
  parents <- vapply(terms, function(term) .or(term$parent, NA_character_),
                    character(1))
  if (!all(parents[-1L] == labels[-length(labels)])) return(NULL)
  groups <- lapply(labels, function(label) {
    droplevels(factor(design$groups[[label]]))
  })
  names(groups) <- labels
  if (any(vapply(groups, function(g) anyNA(g), logical(1)))) return(NULL)
  list(labels = labels, groups = groups, depth = length(labels))
}

# Rows are parents, columns are children, each parent row averaging its own
# children with equal weight.
.equal_weight_map <- function(child_parent, parent_count) {
  M <- matrix(0, parent_count, length(child_parent))
  for (parent in seq_len(parent_count)) {
    children <- which(child_parent == parent)
    if (!length(children)) {
      stop("A nesting level contains an empty unit.", call. = FALSE)
    }
    M[parent, children] <- 1 / length(children)
  }
  M
}

# Composed observation -> level-`target` membership matrix, equal weight at
# every stage of the descent.
.cascade_membership <- function(chain, target) {
  depth <- chain$depth
  deepest <- chain$groups[[depth]]
  M <- .equal_weight_map(as.integer(deepest), nlevels(deepest))
  current <- deepest
  level <- depth - 1L
  while (level >= target) {
    parent <- chain$groups[[level]]
    unit_parent <- vapply(split(as.integer(parent), current), function(values) {
      if (length(unique(values)) != 1L) {
        stop("Every nested unit must belong to exactly one parent.",
             call. = FALSE)
      }
      values[[1L]]
    }, integer(1))
    M <- .equal_weight_map(unit_parent, nlevels(parent)) %*% M
    current <- parent
    level <- level - 1L
  }
  list(M = M, unit = current)
}

# Equal-weight aggregated Gower matrix at level `target`, plus the observation
# index of one representative per unit (for reading unit-level design rows).
.cascade_aggregate <- function(design, chain, target, tol = 1e-10) {
  cascade <- .cascade_membership(chain, target)
  M <- cascade$M
  units <- nrow(M)
  J <- diag(units) - matrix(1 / units, units, units)
  G <- J %*% (M %*% design$G %*% t(M)) %*% J
  G <- (G + t(G)) / 2
  G[abs(G) < tol] <- 0
  unit <- cascade$unit
  blocks <- unname(split(seq_len(design$n), unit))
  list(
    G = G, membership = M, unit = unit, units = units,
    blocks = blocks, sizes = lengths(blocks),
    representatives = vapply(blocks, `[[`, integer(1), 1L)
  )
}

# Recursive exchangeability signature for the units at `level`: a unit's
# signature is the sorted multiset of its children's signatures, bottoming out
# at observation counts. Two units may only be exchanged when their whole
# subtree shape matches, which is the multi-level generalisation of the
# "sorted child observation-count" rule used by the parent/child engine.
.nested_topology_signature <- function(chain, level) {
  depth <- chain$depth
  deepest <- chain$groups[[depth]]
  signature <- as.character(as.integer(table(deepest)))
  current <- deepest
  index <- depth - 1L
  while (index >= level) {
    parent <- chain$groups[[index]]
    unit_parent <- vapply(split(as.integer(parent), current), `[[`, integer(1), 1L)
    signature <- vapply(seq_len(nlevels(parent)), function(p) {
      paste0("(", paste(sort(signature[unit_parent == p]), collapse = ","), ")")
    }, character(1))
    current <- parent
    index <- index - 1L
  }
  signature
}

# Shallowest chain level within which the treatment is constant. Returns NA when
# the treatment varies inside the deepest units (an observation-level term).
.nested_treatment_level <- function(chain, value) {
  constant <- vapply(seq_len(chain$depth), function(level) {
    all(vapply(split(as.character(value), chain$groups[[level]]),
               function(v) length(unique(v)) == 1L, logical(1)))
  }, logical(1))
  if (!any(constant)) return(NA_integer_)
  which(constant)[[1L]]
}


# Assign each fixed term to the shallowest hierarchy level whose units it is
# constant within. A term whose columns are not spanned by any level's unit
# indicators varies inside the deepest unit and belongs to the observation
# stratum; NA marks those for the caller.
.nested_term_levels <- function(design, chain, tol = 1e-10) {
  X <- design$X
  assignment <- attr(X, "assign")
  labels <- attr(stats::terms(design$fixed_formula), "term.labels")
  level_projections <- lapply(chain$groups, function(group) {
    projection_matrix(stats::model.matrix(~ 0 + group), tol = tol)$H
  })
  stats::setNames(vapply(seq_along(labels), function(index) {
    columns <- which(assignment == index)
    if (!length(columns)) return(NA_integer_)
    X_term <- X[, columns, drop = FALSE]
    scale <- max(1, max(abs(X_term)))
    for (level in seq_len(chain$depth)) {
      if (max(abs(level_projections[[level]] %*% X_term - X_term)) <=
          sqrt(tol) * scale) {
        return(level)
      }
    }
    NA_integer_
  }, integer(1)), labels)
}

#' Multi-level unbalanced nested PERMANOVA
#'
#' Tests one or more fixed terms in an unbalanced `/`-generated hierarchy of any
#' depth. The two-level engines
#' ([permanova_unbalanced_nested()] and [permanova_unbalanced_nested_lower()])
#' each test a single categorical treatment in one parent/child pair; this
#' engine covers deeper chains such as `(1 | region/site/quadrat)` and any
#' number of crossed fixed factors, interactions, and covariates.
#'
#' Each fixed term is assigned to the shallowest hierarchy level within which it
#' is constant, and tested on **equally weighted centroids at that level** over
#' the residual among those centroids after the parent level and the other terms
#' sharing the level. Aggregation cascades up the chain -- observations onto the
#' deepest units, those onto their parents with equal weight per child, and so
#' on -- so a level's centroid is not reweighted by how many observations happen
#' to sit beneath it. A term that varies inside the deepest unit (typically an
#' observation-level covariate) is tested over the observation residual with
#' residuals permuted within the deepest units.
#'
#' Permutation moves intact units within their own parent, and only across units
#' whose complete subtree shape matches, so an unequal hierarchy never exchanges
#' structurally different units. No unbalanced expected-mean-square
#' approximation is made.
#'
#' @param design Compiled Design Object with one contiguous `/`-generated
#'   hierarchy of at least two levels.
#' @param permutations A non-negative integer or `"all"`.
#' @param seed Optional random seed.
#' @param hypothesis Sum-of-squares partition applied among the terms sharing a
#'   hierarchy level (`"I"`, `"II"`, or `"III"`).
#' @param keep_permutations Retain the permutation statistics on each test spec.
#' @param max_exact Maximum permutation space allowed for exact enumeration.
#' @param tol Numerical tolerance.
#' @return An object of class `permanova_nested`.
#' @export
permanova_unbalanced_nested_multi <- function(
    design, permutations = 999L, seed = NULL, hypothesis = "I",
    keep_permutations = FALSE, max_exact = 100000L, tol = 1e-10) {
  hypothesis <- match.arg(hypothesis, c("I", "II", "III"))
  chain <- .nested_chain(design)
  if (is.null(chain)) {
    stop("Multi-level unbalanced nested inference requires one contiguous ",
         "`/`-generated hierarchy.", call. = FALSE)
  }
  if (identical(hypothesis, "III")) .validate_type3_contrasts(design, tol = tol)
  X <- design$X
  assignment <- attr(X, "assign")
  labels <- attr(stats::terms(design$fixed_formula), "term.labels")
  if (!length(labels)) {
    stop("At least one fixed-effect term is required.", call. = FALSE)
  }
  term_level <- .nested_term_levels(design, chain, tol = tol)
  if (!is.null(seed)) set.seed(seed)
  tested <- list()
  strata_rows <- list()

  # ---- hierarchy-level terms, level by level -------------------------------
  for (level in seq_len(chain$depth)) {
    level_terms <- labels[which(!is.na(term_level) & term_level == level)]
    if (!length(level_terms)) next
    aggregated <- .cascade_aggregate(design, chain, level, tol = tol)
    units <- aggregated$units
    Xu <- X[aggregated$representatives, , drop = FALSE]
    parent <- if (level > 1L) {
      factor(vapply(split(as.character(chain$groups[[level - 1L]]),
                          aggregated$unit), `[[`, character(1), 1L))
    } else NULL
    parent_columns <- if (is.null(parent)) {
      matrix(1, units, 1L)
    } else {
      stats::model.matrix(~ parent)
    }
    level_columns <- which(assignment %in% match(level_terms, labels))

    stratum_reduced <- projection_matrix(
      cbind(parent_columns, Xu[, level_columns, drop = FALSE]), tol = tol
    )
    stratum_H <- diag(units) - stratum_reduced$H
    stratum_df <- units - stratum_reduced$rank
    if (stratum_df <= 0L) {
      stop("No residual degrees of freedom remain among `",
           chain$labels[[level]], "` units after the terms assigned there.",
           call. = FALSE)
    }
    stratum_ss <- trace_hgh(stratum_H, aggregated$G)
    stratum_ms <- stratum_ss / stratum_df

    signature <- .nested_topology_signature(chain, level)
    strata <- droplevels(interaction(
      if (is.null(parent)) factor(rep("grand mean", units)) else parent,
      signature, drop = TRUE, sep = " | "
    ))
    plan <- .stratified_permutation_plan(strata, permutations, max_exact)
    strata_label <- if (level > 1L) {
      paste0("`", chain$labels[[level - 1L]], "` by identical subtree topology")
    } else "identical subtree topology"

    for (index in seq_along(level_terms)) {
      term <- level_terms[[index]]
      adjust <- level_terms[.ss_type_adjust(index, level_terms, hypothesis)]
      adjust_columns <- which(assignment %in% match(adjust, labels))
      term_columns <- which(assignment == match(term, labels))
      reduced <- projection_matrix(
        cbind(parent_columns, Xu[, adjust_columns, drop = FALSE]), tol = tol
      )
      full <- projection_matrix(
        cbind(parent_columns, Xu[, c(adjust_columns, term_columns), drop = FALSE]),
        tol = tol
      )
      numerator_df <- full$rank - reduced$rank
      if (numerator_df <= 0L) {
        stop("Term `", term, "` is aliased with the hierarchy at level `",
             chain$labels[[level]], "` for Type ", hypothesis, " inference.",
             call. = FALSE)
      }
      numerator_H <- full$H - reduced$H
      contrasts <- .projection_in_subspace(numerator_H, Xu[, term_columns, drop = FALSE],
                                           tol = tol)$projected
      if (.stratified_treatment_rank(contrasts, strata, tol) < numerator_df) {
        stop("Parent-by-topology strata do not span every contrast of `", term,
             "` at level `", chain$labels[[level]],
             "`; intact-unit permutation cannot test the full hypothesis.",
             call. = FALSE)
      }
      numerator_ss <- trace_hgh(numerator_H, aggregated$G)
      numerator_ms <- numerator_ss / numerator_df
      statistic <- numerator_ms / stratum_ms
      statistics <- .freedman_lane_statistics(
        aggregated$G, reduced$H, numerator_H, numerator_df,
        stratum_H, stratum_df, plan$orders,
        label = paste0("Permuting intact `", chain$labels[[level]],
                       "` units for term `", term, "`"), tol = tol
      )
      inference <- .permutation_inference(statistics, statistic, plan$exact, tol)
      p_value <- inference$p_value
      tested[[term]] <- list(
        term = term, kind = "hierarchy", level = chain$labels[[level]],
        numerator_ss = numerator_ss, numerator_df = numerator_df,
        numerator_ms = numerator_ms, numerator_H = numerator_H,
        denominator_name = paste0(chain$labels[[level]], "-centroid Residual"),
        denominator_ms = stratum_ms, denominator_df = stratum_df,
        statistic = statistic, p_value = p_value,
        exchangeable = chain$labels[[level]],
        permutation_strata = strata_label,
        attainable_permutations = plan$attainable,
        distinct_statistics = inference$distinct_statistics,
        minimum_attainable_p = inference$minimum_attainable_p,
        resolution_note = inference$resolution_note,
        permutations_used = nrow(plan$orders),
        exact_permutation = plan$exact && nrow(plan$orders) > 0L,
        permutation_statistics = if (keep_permutations) statistics else NULL
      )
    }
    strata_rows[[length(strata_rows) + 1L]] <- data.frame(
      term = paste0(chain$labels[[level]], " stratum"),
      df = stratum_df, ss = stratum_ss, ms = stratum_ms,
      denominator = NA_character_, denominator_df = NA_real_,
      pseudo_f = NA_real_, p_perm = NA_real_, stringsAsFactors = FALSE
    )
  }

  # ---- observation-level terms over the within-unit residual ---------------
  observation_terms <- labels[is.na(term_level)]
  deepest <- chain$groups[[chain$depth]]
  H_between <- projection_matrix(stats::model.matrix(~ 0 + deepest), tol = tol)$H
  residual_H <- diag(design$n) - H_between
  observation_columns <- which(assignment %in% match(observation_terms, labels))
  within_full <- projection_matrix(
    residual_H %*% X[, observation_columns, drop = FALSE], tol = tol
  )
  residual_df <- (design$n - nlevels(deepest)) - within_full$rank
  residual_H_final <- residual_H - within_full$H
  residual_ss <- trace_hgh(residual_H_final, design$G)
  residual_ms <- if (residual_df > 0L) residual_ss / residual_df else NA_real_

  if (length(observation_terms)) {
    if (residual_df <= 0L) {
      stop("No within-unit residual degrees of freedom remain after the ",
           "observation-level terms.", call. = FALSE)
    }
    for (index in seq_along(observation_terms)) {
      term <- observation_terms[[index]]
      adjust <- observation_terms[
        .ss_type_adjust(index, observation_terms, hypothesis)
      ]
      adjust_columns <- which(assignment %in% match(adjust, labels))
      term_columns <- which(assignment == match(term, labels))
      within_adjust <- residual_H %*% X[, adjust_columns, drop = FALSE]
      within_term <- residual_H %*% X[, term_columns, drop = FALSE]
      reduced <- projection_matrix(within_adjust, tol = tol)
      full <- projection_matrix(cbind(within_adjust, within_term), tol = tol)
      numerator_df <- full$rank - reduced$rank
      if (numerator_df <= 0L) {
        stop("Observation-level term `", term, "` is not estimable for Type ",
             hypothesis, " inference.", call. = FALSE)
      }
      numerator_H <- full$H - reduced$H
      numerator_ss <- trace_hgh(numerator_H, design$G)
      numerator_ms <- numerator_ss / numerator_df
      statistic <- numerator_ms / residual_ms
      result <- .within_unit_freedman_lane(
        design$G, numerator_H, numerator_df, H_between + reduced$H,
        residual_H_final, residual_df, statistic, deepest,
        permutations, max_exact, tol,
        label = paste0("Permuting observations within `",
                       chain$labels[[chain$depth]], "` for term `", term, "`")
      )
      tested[[term]] <- list(
        term = term, kind = "observation",
        level = paste0("within ", chain$labels[[chain$depth]]),
        numerator_ss = numerator_ss, numerator_df = numerator_df,
        numerator_ms = numerator_ms, numerator_H = numerator_H,
        denominator_name = "Residual",
        denominator_ms = residual_ms, denominator_df = residual_df,
        statistic = statistic, p_value = result$p_value,
        exchangeable = paste0("observations within ", chain$labels[[chain$depth]]),
        permutation_strata = chain$labels[[chain$depth]],
        attainable_permutations = result$attainable,
        permutations_used = result$used,
        exact_permutation = result$exact,
        permutation_statistics = if (keep_permutations) result$statistics else NULL
      )
    }
  }

  tested <- tested[intersect(labels, names(tested))]
  specs <- stats::setNames(lapply(tested, function(item) {
    spec <- list(
      term = item$term, kind = item$kind, hierarchy_level = item$level,
      numerator = list(ss = item$numerator_ss, df = item$numerator_df,
                       ms = item$numerator_ms, projector = item$numerator_H),
      denominator = list(expression = paste0("MS[", item$denominator_name, "]"),
                         ms = item$denominator_ms, df = item$denominator_df),
      statistic = item$statistic, p_value = item$p_value,
      estimand = paste0("equally weighted `", item$level,
                        "` centroids from cascaded hierarchy aggregation"),
      exchangeable_units = item$exchangeable,
      permutation_strata = item$permutation_strata,
      attainable_permutations = item$attainable_permutations,
      distinct_statistics = item$distinct_statistics,
      minimum_attainable_p = item$minimum_attainable_p,
      permutations_used = item$permutations_used,
      exact_permutation = item$exact_permutation,
      warnings = item$resolution_note,
      permutation_statistics = item$permutation_statistics
    )
    class(spec) <- "permanova_test_spec"
    spec
  }), names(tested))

  table <- do.call(rbind, lapply(specs, function(spec) data.frame(
    term = spec$term, df = spec$numerator$df, ss = spec$numerator$ss,
    ms = spec$numerator$ms, denominator = spec$denominator$expression,
    denominator_df = spec$denominator$df, pseudo_f = spec$statistic,
    p_perm = spec$p_value, stringsAsFactors = FALSE
  )))
  if (length(strata_rows)) table <- rbind(table, do.call(rbind, strata_rows))
  table <- rbind(table, data.frame(
    term = "Residual", df = residual_df, ss = residual_ss, ms = residual_ms,
    denominator = NA_character_, denominator_df = NA_real_,
    pseudo_f = NA_real_, p_perm = NA_real_, stringsAsFactors = FALSE))
  rownames(table) <- NULL

  out <- list(
    call = match.call(), design = design, table = table, tests = specs,
    hierarchy = chain$labels, term_levels = term_level,
    supported_design = paste(
      "one contiguous unbalanced `/`-generated hierarchy of any depth with any",
      "number of fixed terms; each term is tested on equally weighted cascaded",
      "centroids at its own hierarchy level"
    ),
    permutation = list(
      algorithm = paste("cascaded-centroid Freedman-Lane permutation of intact",
                        "units within parent-by-subtree-topology strata"),
      used = max(vapply(specs, `[[`, numeric(1), "permutations_used"), 0),
      exact = all(vapply(specs, `[[`, logical(1), "exact_permutation")),
      seed = seed
    )
  )
  class(out) <- "permanova_nested"
  out
}
