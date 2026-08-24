.pairwise_contrast_family <- function(structure) {
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
  list(pairs = pairs, contrasts = contrasts)
}

.closed_intersection_subsets <- function(hypotheses, max_intersections) {
  if (length(max_intersections) != 1L || !is.finite(max_intersections) ||
      max_intersections < 1 || max_intersections != as.integer(max_intersections)) {
    stop("`max_intersections` must be a positive integer.", call. = FALSE)
  }
  number <- 2 ^ hypotheses - 1
  if (!is.finite(number) || number > max_intersections) {
    stop(
      "The closed pairwise family requires ", format(number, scientific = FALSE),
      " intersection tests, exceeding `max_intersections = ",
      max_intersections, "`; reduce the number of factor levels or explicitly ",
      "raise the computational guard.",
      call. = FALSE
    )
  }
  unlist(lapply(seq_len(hypotheses), function(size) {
    utils::combn(seq_len(hypotheses), size, simplify = FALSE)
  }), recursive = FALSE)
}

.intersection_components <- function(levels, pairs, subset) {
  parent <- seq_along(levels)
  find_root <- function(index) {
    while (parent[[index]] != index) index <- parent[[index]]
    index
  }
  unite <- function(first, second) {
    first_root <- find_root(first)
    second_root <- find_root(second)
    if (first_root != second_root) parent[[second_root]] <<- first_root
  }
  for (edge in subset) {
    unite(match(pairs[[edge]][[1L]], levels),
          match(pairs[[edge]][[2L]], levels))
  }
  roots <- vapply(seq_along(levels), find_root, integer(1))
  components <- unname(split(levels, roots))
  Filter(function(component) length(component) > 1L, components)
}

.component_assignment_orders <- function(position_sets, n) {
  sizes <- lengths(position_sets)
  pool <- sort(unlist(position_sets, use.names = FALSE))
  total <- .safe_factorial(length(pool)) /
    prod(vapply(sizes, .safe_factorial, numeric(1)))
  results <- vector("list", as.integer(total))
  cursor <- 0L
  walk <- function(index, remaining, order) {
    if (index == length(position_sets)) {
      order[position_sets[[index]]] <- sort(remaining)
      cursor <<- cursor + 1L
      results[[cursor]] <<- order
      return(invisible(NULL))
    }
    choices <- utils::combn(
      remaining, sizes[[index]], simplify = FALSE
    )
    for (choice in choices) {
      next_order <- order
      next_order[position_sets[[index]]] <- sort(choice)
      walk(index + 1L, setdiff(remaining, choice), next_order)
    }
    invisible(NULL)
  }
  walk(1L, pool, seq_len(n))
  do.call(rbind, results)
}

.combine_component_assignments <- function(order_sets, positions, n) {
  result <- matrix(seq_len(n), nrow = 1L)
  for (index in seq_along(order_sets)) {
    local <- order_sets[[index]]
    old_rows <- rep(seq_len(nrow(result)), each = nrow(local))
    local_rows <- rep(seq_len(nrow(local)), times = nrow(result))
    combined <- result[old_rows, , drop = FALSE]
    combined[, positions[[index]]] <-
      local[local_rows, positions[[index]], drop = FALSE]
    result <- combined
  }
  result
}

.intersection_permutation_plan <- function(
    group, components, permutations, max_exact) {
  if (length(max_exact) != 1L || !is.finite(max_exact) || max_exact < 1 ||
      max_exact != as.integer(max_exact)) {
    stop("`max_exact` must be one positive integer.", call. = FALSE)
  }
  n <- length(group)
  position_sets <- lapply(components, function(component) {
    lapply(component, function(level) which(group == level))
  })
  positions <- lapply(position_sets, function(sets) {
    sort(unlist(sets, use.names = FALSE))
  })
  attainable_by_component <- vapply(position_sets, function(sets) {
    sizes <- lengths(sets)
    .safe_factorial(sum(sizes)) /
      prod(vapply(sizes, .safe_factorial, numeric(1)))
  }, numeric(1))
  attainable <- .safe_product(attainable_by_component)
  exact_requested <- is.character(permutations)
  if (exact_requested &&
      (length(permutations) != 1L || permutations != "all")) {
    stop("`permutations` must be a non-negative integer or \"all\".",
         call. = FALSE)
  }
  if (exact_requested && (!is.finite(attainable) || attainable > max_exact)) {
    stop(
      "Exact intersection enumeration exceeds `max_exact`; request a finite Monte Carlo sample.",
      call. = FALSE
    )
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
    order_sets <- lapply(position_sets, .component_assignment_orders, n = n)
    orders <- .combine_component_assignments(order_sets, positions, n)
    return(list(
      orders = orders, exact = TRUE, attainable = attainable,
      attainable_by_component = attainable_by_component,
      positions = positions
    ))
  }
  permutations <- as.integer(permutations)
  orders <- matrix(
    rep(seq_len(n), each = permutations), nrow = permutations,
    ncol = n, byrow = FALSE
  )
  if (permutations) {
    for (row in seq_len(permutations)) {
      for (sets in position_sets) {
        draw <- sample(sort(unlist(sets, use.names = FALSE)))
        offset <- 0L
        for (target in sets) {
          selected <- seq.int(offset + 1L, offset + length(target))
          orders[row, target] <- sort(draw[selected])
          offset <- offset + length(target)
        }
      }
    }
  }
  list(
    orders = orders, exact = FALSE, attainable = attainable,
    attainable_by_component = attainable_by_component,
    positions = positions
  )
}

.strong_pairwise_statistics <- function(
    design, structure, projectors, observed, subset, plan, tol) {
  statistics <- .permuted_statistics(
    nrow(plan$orders),
    function(row) {
      order <- plan$orders[row, ]
      G <- design$G[order, order, drop = FALSE]
      denominator <- trace_hgh(structure$H_denominator, G) /
        structure$denominator_df
      numerators <- vapply(subset, function(index) {
        trace_hgh(projectors[[index]]$H, G)
      }, numeric(1))
      values <- if (.mean_square_is_zero(denominator, G, tol)) {
        ifelse(abs(numerators) <= tol, 0, Inf)
      } else numerators / denominator
      max(values)
    },
    label = "Permuting labels within connected components"
  )
  observed_maximum <- max(observed[subset])
  p_value <- NA_real_
  if (length(statistics)) {
    extreme <- sum(
      statistics >= observed_maximum -
        max(100 * tol, sqrt(.Machine$double.eps)) * max(1, abs(observed_maximum)),
      na.rm = TRUE)
    p_value <- if (plan$exact) extreme / length(statistics) else
      (extreme + 1) / (length(statistics) + 1)
  }
  list(observed = observed_maximum, statistics = statistics, p_value = p_value)
}

#' Strong pairwise inference by intersection-specific closed permutation tests
#'
#' Tests every pairwise hypothesis for a fixed one-factor design. For each
#' intersection null, observations are reassigned only within connected
#' components of the selected pairwise equalities; groups outside those
#' components remain fixed. A maximum pseudo-F tests that intersection, and
#' closure forms each adjusted p-value as the maximum intersection p-value over
#' all intersections containing that pair.
#'
#' The finite-sample strong-FWER proof requires independent observations and
#' equality of the complete observation distributions within every true-null
#' connected component. It is not a proof for equality of centroids alone when
#' dispersions or other distributional features differ. The method does not
#' use or claim subset pivotality.
#'
#' @param design Canonical fixed one-factor Design Object.
#' @param permutations A non-negative integer or `"all"`.
#' @param seed Optional random seed.
#' @param keep_permutations Retain intersection maximum statistics and label
#'   assignment orders.
#' @param max_exact Maximum unique label assignments for any intersection.
#' @param max_intersections Maximum intersection tests allowed by the
#'   computational guard.
#' @param tol Numerical tolerance.
#' @return An object of class `permanova_pairwise_strong`.
#' @export
permanova_pairwise_strong <- function(
    design, permutations = 999L, seed = NULL,
    keep_permutations = FALSE, max_exact = 100000L,
    max_intersections = 1023L, tol = 1e-10) {
  structure <- .validate_contrast_design(design, tol = tol)
  if (length(structure$levels) < 2L) {
    stop("Strong pairwise inference requires at least two factor levels.",
         call. = FALSE)
  }
  family <- .pairwise_contrast_family(structure)
  subsets <- .closed_intersection_subsets(
    length(family$pairs), max_intersections
  )
  projectors <- lapply(seq_len(nrow(family$contrasts)), function(index) {
    .contrast_projector(
      family$contrasts[index, ], structure$group, structure$counts,
      structure$H_full, tol
    )
  })
  denominator <- trace_hgh(structure$H_denominator, design$G) /
    structure$denominator_df
  if (.mean_square_is_zero(denominator, design$G, tol)) {
    stop("The within-group residual mean square is numerically zero.",
         call. = FALSE)
  }
  observed <- vapply(projectors, function(projector) {
    trace_hgh(projector$H, design$G) / denominator
  }, numeric(1))
  names(observed) <- rownames(family$contrasts)
  if (!is.null(seed)) set.seed(seed)
  intersections <- lapply(seq_along(subsets), function(index) {
    subset <- subsets[[index]]
    components <- .intersection_components(
      structure$levels, family$pairs, subset
    )
    plan <- .intersection_permutation_plan(
      structure$group, components, permutations, max_exact
    )
    result <- .strong_pairwise_statistics(
      design, structure, projectors, observed, subset, plan, tol
    )
    list(
      id = index, subset = subset,
      hypotheses = rownames(family$contrasts)[subset],
      components = components, statistic = result$observed,
      p_value = result$p_value, exact = plan$exact,
      attainable = plan$attainable, used = nrow(plan$orders),
      statistics = if (keep_permutations) result$statistics else NULL,
      orders = if (keep_permutations) plan$orders else NULL
    )
  })
  singleton <- vapply(intersections, function(item) {
    if (length(item$subset) == 1L) item$subset else NA_integer_
  }, integer(1))
  raw <- vapply(seq_along(family$pairs), function(index) {
    intersections[[which(singleton == index)]]$p_value
  }, numeric(1))
  adjusted <- vapply(seq_along(family$pairs), function(index) {
    values <- vapply(
      intersections[vapply(intersections, function(item) {
        index %in% item$subset
      }, logical(1))], `[[`, numeric(1), "p_value"
    )
    if (all(is.na(values))) NA_real_ else max(values, na.rm = TRUE)
  }, numeric(1))
  proof <- paste(
    "For every intersection S, the selected equality edges form connected",
    "components. Under equality of the complete observation distributions",
    "within each component, the joint data law is invariant to every",
    "sample-size-preserving raw-label reassignment used for S. Therefore the",
    "non-null groups outside those components remain fixed, and the",
    "intersection maximum-F permutation test is level alpha. Closed testing",
    "rejects an elementary pair only when every containing intersection is",
    "rejected, giving strong FWER control without a subset-pivotality claim."
  )
  notes <- c(
    "Strong control applies to pairwise equality of complete distributions under independent within-component exchangeability, not to centroid equality under heterogeneous distributions.",
    "Every intersection uses its own null-preserving label-permutation group; non-null groups outside its connected components remain fixed."
  )
  if (!isTRUE(design$spectral$euclidean)) {
    notes <- c(notes,
               "Empirical inference retained the complete non-Euclidean Gower spectrum.")
  }
  tests <- lapply(seq_along(family$pairs), function(index) {
    pair <- family$pairs[[index]]
    test <- list(
      term = structure$fixed$name,
      contrast = rownames(family$contrasts)[[index]],
      coefficients = family$contrasts[index, ],
      hypothesis = paste0(
        "the complete observation distributions of `", pair[[1L]],
        "` and `", pair[[2L]], "` are equal"
      ),
      numerator = list(
        ss = observed[[index]] * denominator, df = 1L,
        ms = observed[[index]] * denominator,
        projector = projectors[[index]]$H
      ),
      denominator = list(
        expression = "MS[within-group Residual]",
        ss = denominator * structure$denominator_df,
        df = structure$denominator_df, ms = denominator,
        projector = structure$H_denominator
      ),
      statistic = observed[[index]], p_value = raw[[index]],
      p_adjusted = adjusted[[index]], adjustment = "closed intersection permutation",
      error_control = "strong FWER for the pairwise complete-distribution equality family",
      exchangeable_units = "independent observations",
      permutation_strata = "intersection-specific connected equality components",
      permutation_algorithm = "sample-size-preserving raw-label reassignment within each intersection component",
      warnings = notes
    )
    class(test) <- "permanova_test_spec"
    test
  })
  names(tests) <- rownames(family$contrasts)
  intersection_table <- data.frame(
    intersection = vapply(intersections, function(item) {
      paste(item$hypotheses, collapse = " & ")
    }, character(1)),
    components = vapply(intersections, function(item) {
      paste(vapply(item$components, paste, collapse = " = ", character(1)),
            collapse = "; ")
    }, character(1)),
    pseudo_f_max = vapply(intersections, `[[`, numeric(1), "statistic"),
    p_perm = vapply(intersections, `[[`, numeric(1), "p_value"),
    exact = vapply(intersections, `[[`, logical(1), "exact"),
    attainable = vapply(intersections, `[[`, numeric(1), "attainable"),
    used = vapply(intersections, `[[`, integer(1), "used"),
    stringsAsFactors = FALSE
  )
  table <- data.frame(
    contrast = rownames(family$contrasts), df = 1L,
    ss = observed * denominator,
    denominator_df = structure$denominator_df,
    pseudo_f = observed, p_perm = raw, p_closed = adjusted,
    error_control = "strong FWER for complete-distribution pairwise nulls",
    stringsAsFactors = FALSE
  )
  out <- list(
    call = match.call(), design = design, factor = structure$fixed$name,
    levels = structure$levels, counts = structure$counts,
    contrast_matrix = family$contrasts, table = table, tests = tests,
    intersections = intersection_table,
    intersection_details = if (keep_permutations) intersections else NULL,
    family = "all pairwise complete-distribution equality hypotheses",
    family_size = length(family$pairs),
    adjustment = "closed intersection-specific permutation",
    error_control = "finite-sample strong FWER under the declared complete-distribution exchangeability nulls",
    exchangeability_proof = proof,
    exclusions = paste(
      "centroid-only nulls with unequal distributions, planned non-pairwise",
      "contrasts, mixed or repeated designs, dependent observations, and",
      "intersection families exceeding the computational guard"
    ),
    permutation = list(
      algorithm = "intersection-specific component label reassignment",
      intersections = length(intersections), seed = seed,
      exact = all(intersection_table$exact),
      maximum_attainable = max(intersection_table$attainable),
      total_used = sum(intersection_table$used)
    )
  )
  class(out) <- "permanova_pairwise_strong"
  out
}

#' @export
print.permanova_pairwise_strong <- function(x, ...) {
  cat("Strong closed pairwise distance inference\n")
  cat("Family:", x$family, "\n")
  cat("Error control:", x$error_control, "\n")
  print(x$table, row.names = FALSE)
  invisible(x)
}
