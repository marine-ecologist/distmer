.aggregate_gower_units <- function(design, unit, tol = 1e-10) {
  unit <- droplevels(factor(unit))
  if (length(unit) != design$n || anyNA(unit)) {
    stop("Aggregation units must define every observation.", call. = FALSE)
  }
  blocks <- unname(split(seq_len(design$n), unit))
  M <- matrix(0, nrow = length(blocks), ncol = design$n)
  for (index in seq_along(blocks)) {
    M[index, blocks[[index]]] <- 1 / length(blocks[[index]])
  }
  J <- diag(length(blocks)) - matrix(1 / length(blocks), length(blocks),
                                    length(blocks))
  G <- J %*% (M %*% design$G %*% t(M)) %*% J
  G[abs(G) < tol] <- 0
  list(
    unit = unit, blocks = blocks, sizes = lengths(blocks),
    membership = M, G = G,
    representatives = vapply(blocks, `[[`, integer(1), 1L)
  )
}

.safe_product <- function(x) {
  if (!length(x)) return(1)
  if (any(!is.finite(x)) || sum(log(x)) > log(.Machine$double.xmax)) Inf else prod(x)
}

.label_assignment_space <- function(treatment, strata) {
  spaces <- vapply(unname(split(treatment, strata)), function(values) {
    counts <- as.integer(table(values))
    .safe_factorial(length(values)) / prod(vapply(counts, .safe_factorial,
                                                  numeric(1)))
  }, numeric(1))
  .safe_product(spaces)
}

.combine_stratum_orders <- function(order_sets, positions, n) {
  result <- matrix(seq_len(n), nrow = 1L)
  for (index in seq_along(order_sets)) {
    local <- order_sets[[index]]
    old <- result
    result <- matrix(rep(old, each = nrow(local)), ncol = n, byrow = FALSE)
    local_rows <- local[rep(seq_len(nrow(local)), times = nrow(old)), , drop = FALSE]
    old_rows <- old[rep(seq_len(nrow(old)), each = nrow(local)), , drop = FALSE]
    result <- old_rows
    result[, positions[[index]]] <- matrix(
      positions[[index]][local_rows], nrow = nrow(local_rows)
    )
  }
  result
}

.stratified_permutation_plan <- function(strata, permutations, max_exact) {
  strata <- droplevels(factor(strata))
  positions <- unname(split(seq_along(strata), strata))
  attainable_by_stratum <- vapply(
    lengths(positions), .safe_factorial, numeric(1)
  )
  attainable <- .safe_product(attainable_by_stratum)
  exact_requested <- is.character(permutations)
  if (exact_requested &&
      (length(permutations) != 1L || permutations != "all")) {
    stop("`permutations` must be a non-negative integer or \"all\".",
         call. = FALSE)
  }
  if (exact_requested && (!is.finite(attainable) || attainable > max_exact)) {
    stop("Exact enumeration exceeds `max_exact`; request a finite Monte Carlo sample.",
         call. = FALSE)
  }
  if (!exact_requested &&
      (length(permutations) != 1L || !is.finite(permutations) ||
       permutations < 0 || permutations != as.integer(permutations))) {
    stop("`permutations` must be a non-negative integer or \"all\".",
         call. = FALSE)
  }
  enumerate <- exact_requested ||
    (!exact_requested && permutations > 0L && is.finite(attainable) &&
     permutations >= attainable - 1L && attainable <= max_exact)
  if (enumerate) {
    order_sets <- lapply(lengths(positions), .all_orders)
    orders <- .combine_stratum_orders(order_sets, positions, length(strata))
    return(list(
      orders = orders, exact = TRUE, attainable = attainable,
      strata = positions, attainable_by_stratum = attainable_by_stratum
    ))
  }
  permutations <- as.integer(permutations)
  orders <- if (permutations) {
    matrix(seq_along(strata), nrow = permutations,
           ncol = length(strata), byrow = TRUE)
  } else {
    matrix(integer(), nrow = 0L, ncol = length(strata))
  }
  if (permutations) {
    for (row in seq_len(permutations)) {
      for (position in positions) orders[row, position] <- sample(position)
    }
  }
  list(
    orders = orders, exact = FALSE, attainable = attainable,
    strata = positions, attainable_by_stratum = attainable_by_stratum
  )
}

.unit_freedman_lane <- function(G, H_numerator, numerator_df, H_reduced,
                                H_denominator, denominator_df, plan, tol,
                                label = NULL) {
  numerator_ss <- trace_hgh(H_numerator, G)
  denominator_ss <- trace_hgh(H_denominator, G)
  numerator_ms <- numerator_ss / numerator_df
  denominator_ms <- denominator_ss / denominator_df
  if (!is.finite(denominator_ms) || abs(denominator_ms) <= tol) {
    stop("The unit-level residual mean square is numerically zero.", call. = FALSE)
  }
  observed <- numerator_ms / denominator_ms
  statistics <- .freedman_lane_statistics(
    G, H_reduced, H_numerator, numerator_df, H_denominator, denominator_df,
    plan$orders, label = label, tol = tol
  )
  p_value <- NA_real_
  if (length(statistics)) {
    extreme <- sum(statistics >= observed - tol, na.rm = TRUE)
    p_value <- if (plan$exact) extreme / length(statistics) else
      (extreme + 1) / (length(statistics) + 1)
  }
  list(
    numerator_ss = numerator_ss, numerator_ms = numerator_ms,
    denominator_ss = denominator_ss, denominator_ms = denominator_ms,
    observed = observed, statistics = statistics, p_value = p_value
  )
}

.one_categorical_fixed <- function(design, context) {
  fixed_variables <- all.vars(design$fixed_formula)
  if (length(design$fixed_effects) != 1L || length(fixed_variables) != 1L ||
      !identical(.canonical_term(design$fixed_effects[[1L]]),
                 fixed_variables[[1L]])) {
    stop(context, " supports exactly one categorical fixed treatment factor.",
         call. = FALSE)
  }
  name <- fixed_variables[[1L]]
  value <- design$data[[name]]
  if (!is.factor(value) || nlevels(droplevels(value)) < 2L) {
    stop(context, " requires a categorical treatment with at least two levels.",
         call. = FALSE)
  }
  list(name = name, value = droplevels(value))
}

.validate_unequal_random_intercept <- function(design, tol = 1e-10) {
  if (!inherits(design, "permanova_design") || is.null(design$G)) {
    stop("`design` must be a compiled Design Object containing `distance`.",
         call. = FALSE)
  }
  if (length(design$random_terms) != 1L ||
      design$random_terms[[1L]]$origin != "main") {
    stop("Phase 11 requires exactly one ordinary random-intercept grouping factor.",
         call. = FALSE)
  }
  fixed <- .one_categorical_fixed(design, "Phase 11")
  cluster_name <- design$random_terms[[1L]]$label
  cluster <- droplevels(design$groups[[cluster_name]])
  aggregated <- .aggregate_gower_units(design, cluster, tol = tol)
  if (length(aggregated$blocks) < nlevels(fixed$value) + 1L) {
    stop("Too few random-intercept units retain a residual denominator.",
         call. = FALSE)
  }
  cluster_treatment <- vapply(aggregated$blocks, function(index) {
    values <- unique(as.character(fixed$value[index]))
    if (length(values) != 1L) {
      stop("Treatment must be constant within every random-intercept unit.",
           call. = FALSE)
    }
    values
  }, character(1))
  cluster_treatment <- factor(cluster_treatment, levels = levels(fixed$value))
  X <- stats::model.matrix(~ cluster_treatment)
  full <- projection_matrix(X, tol = tol)
  reduced <- projection_matrix(matrix(1, nrow(X), 1L), tol = tol)
  numerator_df <- full$rank - reduced$rank
  denominator_df <- nrow(X) - full$rank
  if (numerator_df != nlevels(cluster_treatment) - 1L || denominator_df <= 0L) {
    stop("The cluster-level treatment hypothesis is aliased or has no residual degrees of freedom.",
         call. = FALSE)
  }
  size_strata <- factor(aggregated$sizes, levels = sort(unique(aggregated$sizes)))
  exchangeable_rank <- .stratified_treatment_rank(
    X[, -1L, drop = FALSE], size_strata, tol
  )
  if (exchangeable_rank < numerator_df) {
    stop(
      "Cluster-size strata do not span every treatment contrast; unequal-size whole-cluster permutation cannot test the full hypothesis.",
      call. = FALSE
    )
  }
  list(
    fixed = fixed, cluster_name = cluster_name, cluster = cluster,
    aggregated = aggregated, treatment = cluster_treatment,
    size_strata = size_strata,
    attainable_label_assignments = .label_assignment_space(
      cluster_treatment, size_strata
    ),
    X = X,
    H_reduced = reduced$H, H_numerator = full$H - reduced$H,
    numerator_df = numerator_df,
    H_denominator = diag(nrow(X)) - full$H,
    denominator_df = denominator_df
  )
}

#' Unequal-replication random-intercept PERMANOVA
#'
#' Tests one cluster-level categorical treatment using equally weighted random-
#' intercept unit centroids. Freedman--Lane residuals are permuted only among
#' whole clusters having the same observation count. The target is therefore a
#' cluster-weighted location effect; no observation-weighted EMS or parametric
#' F distribution is inferred.
#'
#' @param design Canonical Design Object with one random intercept.
#' @param permutations A non-negative integer or `"all"`.
#' @param seed Optional random seed.
#' @param keep_permutations Retain permutation statistics.
#' @param max_exact Maximum exact permutation space.
#' @param tol Numerical tolerance.
#' @return An object of class `permanova_unequal_random_intercept`.
#' @export
permanova_unequal_random_intercept <- function(
    design, permutations = 999L, seed = NULL,
    keep_permutations = FALSE, max_exact = 100000L, tol = 1e-10) {
  structure <- .validate_unequal_random_intercept(design, tol = tol)
  if (!is.null(seed)) set.seed(seed)
  plan <- .stratified_permutation_plan(
    structure$size_strata, permutations = permutations, max_exact = max_exact
  )
  result <- .unit_freedman_lane(
    structure$aggregated$G, structure$H_numerator,
    structure$numerator_df, structure$H_reduced,
    structure$H_denominator, structure$denominator_df, plan, tol,
    label = "Permuting whole random-intercept clusters"
  )
  notes <- character()
  if (!isTRUE(design$spectral$euclidean)) {
    notes <- c(notes, "Empirical inference retained the complete non-Euclidean Gower spectrum.")
  }
  if (result$numerator_ss < -tol || result$denominator_ss < -tol) {
    notes <- c(notes, "A negative distance-based sum of squares was retained.")
  }
  test <- list(
    term = structure$fixed$name,
    hypothesis = paste0("cluster-weighted centroids are equal across `",
                        structure$fixed$name, "` levels"),
    estimand = "equally weighted random-intercept unit centroids",
    numerator = list(
      ss = result$numerator_ss, df = structure$numerator_df,
      ms = result$numerator_ms, projector = structure$H_numerator
    ),
    denominator = list(
      expression = "MS[cluster-centroid Residual]",
      ss = result$denominator_ss, df = structure$denominator_df,
      ms = result$denominator_ms, projector = structure$H_denominator
    ),
    statistic = result$observed, p_value = result$p_value,
    exchangeable_units = structure$cluster_name,
    permutation_strata = "exact observation count per random-intercept unit",
    permutation_algorithm = paste(
      "intercept-only cluster-centroid Freedman-Lane residual permutation",
      "within equal-size strata (algebraically identical to raw centroid permutation)"
    ),
    attainable_permutations = plan$attainable,
    attainable_label_assignments = structure$attainable_label_assignments,
    permutations_used = nrow(plan$orders),
    unique_permutations_used = if (nrow(plan$orders)) nrow(unique(plan$orders)) else 0L,
    exact_permutation = plan$exact && nrow(plan$orders) > 0L,
    approximation_status = if (!nrow(plan$orders)) "not tested" else if (plan$exact) {
      "exhaustive size-stratified whole-cluster residual permutation"
    } else "Monte Carlo size-stratified whole-cluster residual permutation",
    warnings = notes,
    permutation_statistics = if (keep_permutations) result$statistics else NULL
  )
  class(test) <- "permanova_test_spec"
  table <- data.frame(
    term = c(structure$fixed$name, "Cluster-centroid Residual"),
    df = c(structure$numerator_df, structure$denominator_df),
    ss = c(result$numerator_ss, result$denominator_ss),
    ms = c(result$numerator_ms, result$denominator_ms),
    pseudo_f = c(result$observed, NA_real_),
    p_perm = c(result$p_value, NA_real_)
  )
  out <- list(
    call = match.call(), design = design, table = table, test = test,
    structure = list(
      treatment = structure$fixed$name,
      random_intercept = structure$cluster_name,
      units = length(structure$aggregated$blocks),
      unit_sizes = structure$aggregated$sizes,
      size_strata = table(structure$size_strata)
    ),
    estimand = test$estimand,
    supported_design = paste(
      "one cluster-level categorical treatment and one random intercept;",
      "cluster-weighted centroid inference conditional on exact cluster size"
    ),
    permutation = list(
      algorithm = test$permutation_algorithm,
      exact = test$exact_permutation, attainable = plan$attainable,
      attainable_label_assignments = structure$attainable_label_assignments,
      used = nrow(plan$orders), seed = seed
    )
  )
  class(out) <- "permanova_unequal_random_intercept"
  out
}

#' @export
print.permanova_unequal_random_intercept <- function(x, ...) {
  cat("Unequal-replication random-intercept PERMANOVA\n")
  cat("Estimand:", x$estimand, "\n")
  cat("Permutation:", x$permutation$algorithm, "\n")
  print(x$table, row.names = FALSE)
  invisible(x)
}

.stratified_treatment_rank <- function(X_treatment, strata, tol) {
  X_treatment <- as.matrix(X_treatment)
  if (!nrow(X_treatment) || !ncol(X_treatment) ||
      length(strata) != nrow(X_treatment)) {
    stop("Treatment columns and exchangeability strata do not conform.",
         call. = FALSE)
  }
  center <- matrix(0, nrow(X_treatment), nrow(X_treatment))
  for (indices in unname(split(seq_len(nrow(X_treatment)), strata))) {
    center[indices, indices] <- diag(length(indices)) -
      matrix(1 / length(indices), length(indices), length(indices))
  }
  projected <- center %*% X_treatment
  original_singular <- svd(X_treatment, nu = 0L, nv = 0L)$d
  projected_singular <- svd(projected, nu = 0L, nv = 0L)$d
  original_scale <- max(original_singular)
  relative_floor <- max(100 * tol, sqrt(.Machine$double.eps))
  absolute_floor <- relative_floor * original_scale
  rank <- sum(projected_singular > absolute_floor)
  attr(rank, "projected_singular_values") <- projected_singular
  attr(rank, "absolute_floor") <- absolute_floor
  attr(rank, "relative_floor") <- relative_floor
  rank
}

.validate_unbalanced_nested <- function(design, tol = 1e-10) {
  if (!inherits(design, "permanova_design") || is.null(design$G)) {
    stop("`design` must be a compiled Design Object containing `distance`.",
         call. = FALSE)
  }
  fixed <- .one_categorical_fixed(design, "Unbalanced nested inference")
  if (length(design$random_terms) != 2L) {
    stop("Unbalanced nested inference requires exactly one parent/child hierarchy.",
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
  aggregated <- .aggregate_gower_units(design, parent, tol = tol)
  parent_treatment <- vapply(aggregated$blocks, function(index) {
    values <- unique(as.character(fixed$value[index]))
    if (length(values) != 1L) {
      stop("Treatment must be constant within every top-level hierarchy unit.",
           call. = FALSE)
    }
    values
  }, character(1))
  parent_treatment <- factor(parent_treatment, levels = levels(fixed$value))

  child_mapping <- unique(data.frame(
    child = as.integer(child), parent = as.integer(parent)
  ))
  if (any(table(child_mapping$child) != 1L)) {
    stop("Every nested child must belong to exactly one parent.", call. = FALSE)
  }
  child_sizes <- as.integer(table(child))
  topology <- vapply(seq_along(aggregated$blocks), function(parent_index) {
    children <- child_mapping$child[child_mapping$parent == parent_index]
    paste(sort(child_sizes[children]), collapse = ",")
  }, character(1))
  topology_strata <- factor(topology, levels = unique(topology))

  X <- stats::model.matrix(~ parent_treatment)
  full <- projection_matrix(X, tol = tol)
  reduced <- projection_matrix(matrix(1, nrow(X), 1L), tol = tol)
  numerator_df <- full$rank - reduced$rank
  denominator_df <- nrow(X) - full$rank
  if (numerator_df != nlevels(parent_treatment) - 1L || denominator_df <= 0L) {
    stop("The top-level treatment hypothesis is aliased or lacks residual degrees of freedom.",
         call. = FALSE)
  }
  if (.stratified_treatment_rank(X[, -1L, drop = FALSE], topology_strata,
                                 tol) < numerator_df) {
    stop(
      "Hierarchy-topology strata do not span every treatment contrast; intact top-unit permutation cannot test the full hypothesis.",
      call. = FALSE
    )
  }
  list(
    fixed = fixed, parent_name = parent_name, child_name = child_name,
    parent = parent, child = child, aggregated = aggregated,
    treatment = parent_treatment, topology = topology,
    topology_strata = topology_strata,
    H_reduced = reduced$H, H_numerator = full$H - reduced$H,
    numerator_df = numerator_df,
    H_denominator = diag(nrow(X)) - full$H,
    denominator_df = denominator_df
  )
}

#' Unbalanced nested top-unit PERMANOVA
#'
#' Tests one treatment assigned to top-level units of one unbalanced nested
#' hierarchy. Equally weighted top-unit centroids are permuted only among
#' hierarchies with the same sorted vector of child observation counts.
#'
#' @inheritParams permanova_unequal_random_intercept
#' @return An object of class `permanova_unbalanced_nested`.
#' @export
permanova_unbalanced_nested <- function(
    design, permutations = 999L, seed = NULL,
    keep_permutations = FALSE, max_exact = 100000L, tol = 1e-10) {
  structure <- .validate_unbalanced_nested(design, tol = tol)
  if (!is.null(seed)) set.seed(seed)
  plan <- .stratified_permutation_plan(
    structure$topology_strata, permutations, max_exact
  )
  result <- .unit_freedman_lane(
    structure$aggregated$G, structure$H_numerator,
    structure$numerator_df, structure$H_reduced,
    structure$H_denominator, structure$denominator_df, plan, tol,
    label = "Permuting intact top-level hierarchy units"
  )
  warnings <- if (!isTRUE(design$spectral$euclidean)) {
    "Empirical inference retained the complete non-Euclidean Gower spectrum."
  } else character()
  test <- list(
    term = structure$fixed$name,
    hypothesis = paste0("top-unit-weighted centroids are equal across `",
                        structure$fixed$name, "` levels"),
    estimand = "equally weighted top-level hierarchy centroids",
    numerator = list(
      ss = result$numerator_ss, df = structure$numerator_df,
      ms = result$numerator_ms, projector = structure$H_numerator
    ),
    denominator = list(
      expression = "MS[top-unit-centroid Residual]",
      ss = result$denominator_ss, df = structure$denominator_df,
      ms = result$denominator_ms, projector = structure$H_denominator
    ),
    statistic = result$observed, p_value = result$p_value,
    exchangeable_units = structure$parent_name,
    descendant_units = c(structure$child_name, "observations"),
    permutation_strata = "identical sorted child observation-count topology",
    permutation_algorithm = "top-centroid Freedman-Lane residual permutation within identical hierarchy topologies",
    attainable_permutations = plan$attainable,
    permutations_used = nrow(plan$orders),
    unique_permutations_used = if (nrow(plan$orders)) nrow(unique(plan$orders)) else 0L,
    exact_permutation = plan$exact && nrow(plan$orders) > 0L,
    approximation_status = if (!nrow(plan$orders)) "not tested" else if (plan$exact) {
      "exhaustive topology-stratified top-unit residual permutation"
    } else "Monte Carlo topology-stratified top-unit residual permutation",
    warnings = warnings,
    permutation_statistics = if (keep_permutations) result$statistics else NULL
  )
  class(test) <- "permanova_test_spec"
  table <- data.frame(
    term = c(structure$fixed$name, "Top-unit-centroid Residual"),
    df = c(structure$numerator_df, structure$denominator_df),
    ss = c(result$numerator_ss, result$denominator_ss),
    ms = c(result$numerator_ms, result$denominator_ms),
    pseudo_f = c(result$observed, NA_real_),
    p_perm = c(result$p_value, NA_real_)
  )
  out <- list(
    call = match.call(), design = design, table = table, test = test,
    structure = list(
      treatment = structure$fixed$name, parent = structure$parent_name,
      child = structure$child_name,
      top_units = length(structure$aggregated$blocks),
      topology_signatures = structure$topology
    ),
    estimand = test$estimand,
    supported_design = paste(
      "one top-level treatment and one unbalanced parent/child hierarchy;",
      "top-unit-weighted inference conditional on complete hierarchy topology"
    ),
    permutation = list(
      algorithm = test$permutation_algorithm,
      exact = test$exact_permutation, attainable = plan$attainable,
      used = nrow(plan$orders), seed = seed
    )
  )
  class(out) <- "permanova_unbalanced_nested"
  out
}

#' @export
print.permanova_unbalanced_nested <- function(x, ...) {
  cat("Unbalanced nested top-unit PERMANOVA\n")
  cat("Estimand:", x$estimand, "\n")
  cat("Permutation:", x$permutation$algorithm, "\n")
  print(x$table, row.names = FALSE)
  invisible(x)
}

.validate_unbalanced_block <- function(design, tol = 1e-10) {
  if (!inherits(design, "permanova_design") || is.null(design$G)) {
    stop("`design` must be a compiled Design Object containing `distance`.",
         call. = FALSE)
  }
  fixed <- .one_categorical_fixed(design, "Unbalanced block inference")
  if (length(design$random_terms) != 2L) {
    stop(
      "Unbalanced block inference requires one block random intercept and one explicit treatment:block interaction.",
      call. = FALSE
    )
  }
  origins <- vapply(design$random_terms, `[[`, character(1), "origin")
  block_terms <- design$random_terms[origins == "main"]
  interaction_terms <- design$random_terms[origins == "interaction"]
  if (length(block_terms) != 1L || length(interaction_terms) != 1L ||
      length(block_terms[[1L]]$variables) != 1L ||
      !setequal(interaction_terms[[1L]]$variables,
                c(fixed$name, block_terms[[1L]]$variables[[1L]]))) {
    stop(
      "Declare exactly `(1 | block)` and explicit `(1 | treatment:block)` random terms.",
      call. = FALSE
    )
  }
  block_name <- block_terms[[1L]]$label
  cell_name <- interaction_terms[[1L]]$label
  block <- droplevels(design$groups[[block_name]])
  cell <- droplevels(design$groups[[cell_name]])
  aggregated <- .aggregate_gower_units(design, cell, tol = tol)
  cell_mapping <- do.call(rbind, lapply(seq_along(aggregated$blocks), function(i) {
    index <- aggregated$blocks[[i]]
    treatment_values <- unique(as.character(fixed$value[index]))
    block_values <- unique(as.character(block[index]))
    if (length(treatment_values) != 1L || length(block_values) != 1L) {
      stop("Every treatment:block cell must map to one treatment and one block.",
           call. = FALSE)
    }
    data.frame(
      cell = i, treatment = treatment_values, block = block_values,
      size = length(index), stringsAsFactors = FALSE
    )
  }))
  if (anyDuplicated(cell_mapping[c("treatment", "block")])) {
    stop("Every observed treatment:block combination must form one cell.",
         call. = FALSE)
  }
  cell_treatment <- factor(cell_mapping$treatment, levels = levels(fixed$value))
  cell_block <- factor(cell_mapping$block, levels = levels(block))
  full_X <- stats::model.matrix(~ cell_treatment + cell_block)
  reduced_X <- stats::model.matrix(~ cell_block)
  full <- projection_matrix(full_X, tol = tol)
  reduced <- projection_matrix(reduced_X, tol = tol)
  numerator_df <- full$rank - reduced$rank
  denominator_df <- nrow(full_X) - full$rank
  if (numerator_df != nlevels(cell_treatment) - 1L || denominator_df <= 0L) {
    stop(
      "The incomplete block incidence is disconnected, aliases a treatment contrast, or leaves no treatment:block residual df.",
      call. = FALSE
    )
  }
  exchangeability_strata <- interaction(
    cell_block, factor(cell_mapping$size), drop = TRUE, lex.order = TRUE
  )
  treatment_columns <- stats::model.matrix(~ 0 + cell_treatment)
  if (.stratified_treatment_rank(
        treatment_columns, exchangeability_strata, tol
      ) < numerator_df) {
    stop(
      "Within-block cell-size strata do not span every treatment contrast; intact cell-profile permutation cannot test the full hypothesis.",
      call. = FALSE
    )
  }
  list(
    fixed = fixed, block_name = block_name, cell_name = cell_name,
    block = block, cell = cell, aggregated = aggregated,
    cell_mapping = cell_mapping, treatment = cell_treatment,
    block_factor = cell_block, exchangeability_strata = exchangeability_strata,
    H_reduced = reduced$H, H_numerator = full$H - reduced$H,
    numerator_df = numerator_df,
    H_denominator = diag(nrow(full_X)) - full$H,
    denominator_df = denominator_df
  )
}

.unit_raw_permutation <- function(G, H_numerator, numerator_df,
                                  H_denominator, denominator_df,
                                  plan, tol) {
  numerator_ss <- trace_hgh(H_numerator, G)
  denominator_ss <- trace_hgh(H_denominator, G)
  numerator_ms <- numerator_ss / numerator_df
  denominator_ms <- denominator_ss / denominator_df
  if (!is.finite(denominator_ms) || abs(denominator_ms) <= tol) {
    stop("The cell-profile residual mean square is numerically zero.",
         call. = FALSE)
  }
  observed <- numerator_ms / denominator_ms
  statistics <- vapply(seq_len(nrow(plan$orders)), function(row) {
    index <- plan$orders[row, ]
    G_permuted <- G[index, index, drop = FALSE]
    (trace_hgh(H_numerator, G_permuted) / numerator_df) /
      (trace_hgh(H_denominator, G_permuted) / denominator_df)
  }, numeric(1))
  p_value <- NA_real_
  if (length(statistics)) {
    extreme <- sum(statistics >= observed - tol, na.rm = TRUE)
    p_value <- if (plan$exact) extreme / length(statistics) else
      (extreme + 1) / (length(statistics) + 1)
  }
  list(
    numerator_ss = numerator_ss, numerator_ms = numerator_ms,
    denominator_ss = denominator_ss, denominator_ms = denominator_ms,
    observed = observed, statistics = statistics, p_value = p_value
  )
}

#' Unbalanced randomized-block cell-profile PERMANOVA
#'
#' Tests one treatment in a randomized block design with incomplete treatment
#' incidence and/or unequal replicated cell sizes. Treatment:block cell
#' centroids receive equal weight and move only within block and exact cell-size
#' strata. The residual treatment:block mean square is descriptive; inference
#' is empirical.
#'
#' @inheritParams permanova_unequal_random_intercept
#' @return An object of class `permanova_unbalanced_block`.
#' @export
permanova_unbalanced_block <- function(
    design, permutations = 999L, seed = NULL,
    keep_permutations = FALSE, max_exact = 100000L, tol = 1e-10) {
  structure <- .validate_unbalanced_block(design, tol = tol)
  if (!is.null(seed)) set.seed(seed)
  plan <- .stratified_permutation_plan(
    structure$exchangeability_strata, permutations, max_exact
  )
  result <- .unit_raw_permutation(
    structure$aggregated$G, structure$H_numerator,
    structure$numerator_df, structure$H_denominator,
    structure$denominator_df, plan, tol
  )
  warnings <- if (!isTRUE(design$spectral$euclidean)) {
    "Empirical inference retained the complete non-Euclidean Gower spectrum."
  } else character()
  test <- list(
    term = structure$fixed$name,
    hypothesis = paste0("cell-profile-weighted centroids are equal across `",
                        structure$fixed$name, "` levels after block adjustment"),
    estimand = "equally weighted observed treatment-by-block cell centroids",
    numerator = list(
      ss = result$numerator_ss, df = structure$numerator_df,
      ms = result$numerator_ms, projector = structure$H_numerator
    ),
    denominator = list(
      expression = "MS[cell-centroid treatment:block Residual]",
      ss = result$denominator_ss, df = structure$denominator_df,
      ms = result$denominator_ms, projector = structure$H_denominator
    ),
    statistic = result$observed, p_value = result$p_value,
    exchangeable_units = structure$cell_name,
    permutation_strata = paste0(structure$block_name,
                                " crossed with exact cell observation count"),
    permutation_algorithm = "raw permutation of intact treatment:block cell centroids within block-by-size strata",
    attainable_permutations = plan$attainable,
    permutations_used = nrow(plan$orders),
    unique_permutations_used = if (nrow(plan$orders)) nrow(unique(plan$orders)) else 0L,
    exact_permutation = plan$exact && nrow(plan$orders) > 0L,
    approximation_status = if (!nrow(plan$orders)) "not tested" else if (plan$exact) {
      "exhaustive within-block-by-size cell-profile permutation"
    } else "Monte Carlo within-block-by-size cell-profile permutation",
    warnings = warnings,
    permutation_statistics = if (keep_permutations) result$statistics else NULL
  )
  class(test) <- "permanova_test_spec"
  table <- data.frame(
    term = c(structure$fixed$name, "Cell-centroid treatment:block Residual"),
    df = c(structure$numerator_df, structure$denominator_df),
    ss = c(result$numerator_ss, result$denominator_ss),
    ms = c(result$numerator_ms, result$denominator_ms),
    pseudo_f = c(result$observed, NA_real_),
    p_perm = c(result$p_value, NA_real_)
  )
  out <- list(
    call = match.call(), design = design, table = table, test = test,
    structure = list(
      treatment = structure$fixed$name, block = structure$block_name,
      cell = structure$cell_name, cells = nrow(structure$cell_mapping),
      observed_incidence = table(
        structure$cell_mapping$block, structure$cell_mapping$treatment
      ),
      cell_sizes = structure$cell_mapping$size
    ),
    estimand = test$estimand,
    supported_design = paste(
      "one treatment crossed with one block factor and an explicit",
      "treatment:block cell component; cell-weighted empirical inference"
    ),
    permutation = list(
      algorithm = test$permutation_algorithm,
      exact = test$exact_permutation, attainable = plan$attainable,
      used = nrow(plan$orders), seed = seed
    )
  )
  class(out) <- "permanova_unbalanced_block"
  out
}

#' @export
print.permanova_unbalanced_block <- function(x, ...) {
  cat("Unbalanced randomized-block cell-profile PERMANOVA\n")
  cat("Estimand:", x$estimand, "\n")
  cat("Permutation:", x$permutation$algorithm, "\n")
  print(x$table, row.names = FALSE)
  invisible(x)
}
