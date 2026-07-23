.spatial_median <- function(points, tol = 1e-10) {
  points <- as.matrix(points)
  if (!ncol(points)) return(numeric())
  if (nrow(points) == 1L) return(as.numeric(points[1L, ]))
  objective <- function(candidate) {
    sum(sqrt(rowSums(sweep(points, 2L, candidate, "-") ^ 2)))
  }
  gradient <- function(candidate) {
    delta <- sweep(points, 2L, candidate, "-")
    lengths <- sqrt(rowSums(delta ^ 2))
    nonzero <- lengths > tol
    if (!any(nonzero)) return(rep(0, ncol(points)))
    -colSums(delta[nonzero, , drop = FALSE] / lengths[nonzero])
  }
  start <- apply(points, 2L, stats::median)
  as.numeric(stats::optim(start, objective, gradient, method = "BFGS")$par)
}

.dispersion_geometry <- function(distance, group,
                                 type = c("median", "centroid"),
                                 bias_adjust = FALSE, tol = 1e-10) {
  type <- match.arg(type)
  D <- as_distance_matrix(distance)
  group <- droplevels(factor(group))
  if (length(group) != nrow(D)) {
    stop("Dispersion groups must match the distance-matrix observations.",
         call. = FALSE)
  }
  if (nlevels(group) < 2L) {
    stop("Dispersion inference requires at least two groups.", call. = FALSE)
  }
  group_sizes <- table(group)
  if (any(group_sizes < 2L)) {
    stop("Every dispersion group must contain at least two observations.",
         call. = FALSE)
  }

  G <- gower_center(D)
  decomposition <- eigen(G, symmetric = TRUE)
  eig <- decomposition$values
  spectral_tol <- sqrt(.Machine$double.eps)
  threshold <- max(spectral_tol, spectral_tol * max(abs(eig)))
  keep <- abs(eig) > threshold
  eig_retained <- eig[keep]
  coordinates <- decomposition$vectors[, keep, drop = FALSE]
  if (length(eig_retained)) {
    coordinates <- coordinates %*% diag(
      sqrt(abs(eig_retained)), nrow = length(eig_retained)
    )
  }
  positive <- eig_retained > 0
  groups <- levels(group)
  centers <- matrix(
    0, nrow = length(groups), ncol = ncol(coordinates),
    dimnames = list(groups, paste0("PCoA", seq_len(ncol(coordinates))))
  )
  center_block <- function(columns) {
    if (!length(columns)) return(invisible(NULL))
    for (level in groups) {
      block <- coordinates[group == level, columns, drop = FALSE]
      centers[level, columns] <<- if (type == "centroid") {
        colMeans(block)
      } else {
        .spatial_median(block, tol = tol)
      }
    }
    invisible(NULL)
  }
  center_block(which(positive))
  center_block(which(!positive))

  squared_from_center <- function(columns) {
    if (!length(columns)) return(rep(0, nrow(D)))
    delta <- coordinates[, columns, drop = FALSE] -
      centers[as.character(group), columns, drop = FALSE]
    rowSums(delta ^ 2)
  }
  positive_squared <- squared_from_center(which(positive))
  negative_squared <- squared_from_center(which(!positive))
  signed_squared <- positive_squared - negative_squared
  clamped <- signed_squared < -tol
  near_zero <- signed_squared < 0 & !clamped
  signed_squared[clamped | near_zero] <- 0
  distances <- sqrt(signed_squared)
  if (isTRUE(bias_adjust)) {
    distances <- distances * sqrt(
      as.numeric(group_sizes[group]) / (as.numeric(group_sizes[group]) - 1)
    )
  }
  names(distances) <- rownames(D)
  group_means <- tapply(distances, group, mean)

  list(
    distances = distances, group = group, group_distances = group_means,
    centers = centers, coordinates = coordinates,
    eigenvalues = eig,
    retained_eigenvalues = setNames(
      eig_retained, paste0("PCoA", seq_along(eig_retained))
    ),
    positive_axes = sum(positive), negative_axes = sum(!positive),
    negative_count = sum(eig < -tol),
    zero_clamped = sum(clamped),
    type = type, bias_adjust = isTRUE(bias_adjust),
    warnings = if (any(clamped)) {
      paste0(sum(clamped),
             " negative squared distance(s) to the group center were changed to zero.")
    } else character()
  )
}

.validate_dispersion_target <- function(design, term) {
  if (!inherits(design, "permanova_design")) {
    stop("`design` must be produced by `compile_permanova_design()`.",
         call. = FALSE)
  }
  if (is.null(design$distance) || is.null(design$G)) {
    stop("The Design Object must include `distance`.", call. = FALSE)
  }
  if (length(term) != 1L || !is.character(term) || is.na(term) || !nzchar(term)) {
    stop("`term` must name one categorical fixed term.", call. = FALSE)
  }
  if (length(design$fixed_effects) != 1L ||
      !identical(term, design$fixed_effects[[1L]])) {
    stop(
      "Phase 8 dispersion inference requires a design containing exactly the one requested fixed term.",
      call. = FALSE
    )
  }
  variables <- all.vars(design$fixed_formula)
  if (length(variables) != 1L || !identical(.canonical_term(term), variables[[1L]]) ||
      !is.factor(design$data[[variables[[1L]]]])) {
    stop("The Phase 8 dispersion target must be one categorical main-effect factor.",
         call. = FALSE)
  }
  group <- droplevels(design$data[[variables[[1L]]]])
  sizes <- table(group)
  if (nlevels(group) < 2L || any(sizes < 2L)) {
    stop("Every dispersion group must contain at least two observations.",
         call. = FALSE)
  }
  if (length(unique(as.integer(sizes))) != 1L) {
    stop("Phase 8 dispersion inference requires equal observations per target group.",
         call. = FALSE)
  }
  list(term = term, variable = variables[[1L]], group = group,
       group_sizes = sizes)
}

.dispersion_response_design <- function(design, distances, tol) {
  response <- design
  response$distance <- as_distance_matrix(stats::dist(distances))
  response$G <- gower_center(response$distance)
  response$spectral <- distance_diagnostics(response$distance, tol = tol)
  response$response <- "distance_to_group_center"
  response$status <- "dispersion-distance response compiled"
  response
}

.dispersion_plan <- function(design, response_design, permutations,
                             max_exact, tol) {
  origins <- if (length(design$random_terms)) {
    vapply(design$random_terms, `[[`, character(1), "origin")
  } else character()
  term <- design$fixed_effects[[1L]]
  projection <- response_design$projections$terms[[term]]

  if (!length(origins)) {
    if (response_design$estimability$fixed_rank >= response_design$n) {
      stop("No residual degrees of freedom remain for dispersion inference.",
           call. = FALSE)
    }
    plan <- .permutation_orders(response_design$n, permutations, max_exact)
    unit_blocks <- lapply(seq_len(response_design$n), function(i) i)
    return(list(
      engine = "one-way", plan = plan, unit_blocks = unit_blocks,
      numerator = list(H = projection$H_term, df = projection$df,
                       ems = c(residual = 1)),
      denominator = list(
        H = diag(response_design$n) - response_design$projections$fixed$H,
        df = response_design$n - response_design$projections$fixed$rank,
        term = "Residual", ems = c(residual = 1)
      ),
      exchangeable_units = "observations",
      permutation_strata = "unrestricted",
      algorithm = "full-model residual permutation of observations"
    ))
  }

  if (length(origins) == 1L && identical(origins, "main")) {
    structure <- .validate_balanced_random_intercept(response_design, tol = tol)
    plan <- .permutation_orders(structure$groups, permutations, max_exact)
    ems <- derive_ems(response_design, tol = tol)[[term]]
    return(list(
      engine = "balanced random intercept", plan = plan,
      unit_blocks = structure$blocks,
      numerator = list(H = projection$H_term, df = projection$df,
                       ems = ems$numerator_under_null),
      denominator = list(H = structure$H_denominator,
                         df = structure$denominator_df,
                         term = paste0(structure$cluster_name,
                                       " within fixed model"),
                         ems = ems$denominator),
      exchangeable_units = structure$cluster_name,
      permutation_strata = "unrestricted whole plots",
      algorithm = "full-model residual permutation of whole plots"
    ))
  }

  if (length(origins) >= 2L && identical(origins[[1L]], "main") &&
      all(origins[-1L] == "nesting")) {
    structure <- .validate_nested_design(response_design, tol = tol)
    level <- structure$term_levels[[term]]
    plan <- .nested_permutation_plan(
      structure, level, permutations = permutations, max_exact = max_exact
    )
    ems <- derive_nested_ems(response_design, tol = tol)[[term]]
    stratum <- structure$strata[[level]]
    return(list(
      engine = "balanced nested hierarchy", plan = plan,
      unit_blocks = plan$unit_blocks,
      numerator = list(H = projection$H_term, df = projection$df,
                       ems = ems$numerator_under_null),
      denominator = list(H = stratum$H, df = stratum$df,
                         term = stratum$term, ems = ems$denominator),
      exchangeable_units = plan$level_name,
      permutation_strata = plan$parent_name,
      algorithm = paste0(
        "full-model residual permutation of `", plan$level_name,
        "` units within `", plan$parent_name,
        "`, moving descendants intact"
      )
    ))
  }

  if (length(origins) == 2L && setequal(origins, c("main", "interaction"))) {
    structure <- .validate_crossed_design(
      response_design, assume_additive = FALSE, tol = tol
    )
    permutation_structure <- .validate_crossed_permutation_design(
      response_design, structure, assume_additive = FALSE, tol = tol
    )
    plan <- .crossed_block_permutation_plan(
      permutation_structure, permutations = permutations,
      max_exact = max_exact
    )
    return(list(
      engine = "balanced crossed randomized block", plan = plan,
      unit_blocks = permutation_structure$cell_blocks,
      numerator = permutation_structure$numerator,
      denominator = c(permutation_structure$denominator,
                       list(term = permutation_structure$denominator$term)),
      exchangeable_units = permutation_structure$interaction_name,
      permutation_strata = permutation_structure$block_name,
      algorithm = paste0(
        "full-model residual permutation of intact `",
        permutation_structure$interaction_name, "` cells within `",
        permutation_structure$block_name, "`"
      )
    ))
  }

  stop(
    "Dispersion inference does not support this random-effect structure; use a validated balanced one-way, single-intercept, nested, or canonical randomized-block design.",
    call. = FALSE
  )
}

.dispersion_statistics <- function(response, response_design, setup, tol) {
  numerator <- setup$numerator
  denominator <- setup$denominator
  numerator_ss <- trace_hgh(numerator$H, response_design$G)
  denominator_ss <- trace_hgh(denominator$H, response_design$G)
  numerator_ms <- numerator_ss / numerator$df
  denominator_ms <- denominator_ss / denominator$df
  if (!is.finite(denominator_ms) || abs(denominator_ms) <= tol) {
    stop("The dispersion-test denominator mean square is numerically zero.",
         call. = FALSE)
  }
  observed <- numerator_ms / denominator_ms

  full_residuals <- as.vector(
    (diag(response_design$n) - response_design$projections$fixed$H) %*% response
  )
  permuted <- numeric(nrow(setup$plan$orders))
  for (row in seq_len(nrow(setup$plan$orders))) {
    index <- .block_permutation_index(
      setup$unit_blocks, setup$plan$orders[row, ]
    )
    permuted_response <- full_residuals[index]
    G_permuted <- tcrossprod(permuted_response)
    numerator_perm <- trace_hgh(numerator$H, G_permuted) / numerator$df
    denominator_perm <- trace_hgh(denominator$H, G_permuted) / denominator$df
    permuted[[row]] <- numerator_perm / denominator_perm
  }
  p_value <- NA_real_
  if (length(permuted)) {
    extreme <- sum(permuted >= observed - tol, na.rm = TRUE)
    p_value <- if (setup$plan$exact) {
      extreme / length(permuted)
    } else {
      (extreme + 1) / (length(permuted) + 1)
    }
  }
  list(
    observed = observed, p_value = p_value, statistics = permuted,
    numerator_ss = numerator_ss, numerator_ms = numerator_ms,
    denominator_ss = denominator_ss, denominator_ms = denominator_ms
  )
}

#' Design-aware multivariate dispersion inference
#'
#' Converts the original dissimilarity to distances from each observation to
#' its target-group spatial median or centroid, then tests that scalar response
#' with the EMS denominator and restricted residual permutation appropriate to
#' a validated balanced design. The initial envelope contains exactly one
#' categorical fixed main effect.
#'
#' @param formula A model formula, or a compiled `permanova_design` object. The
#'   dissimilarity may be supplied on the left-hand side, for example
#'   `bray ~ treatment + (1 | site)`, or through `distance` with a one-sided
#'   formula. A compiled design remains accepted for low-level workflows.
#' @param data A data frame with one row per observation. Do not supply it when
#'   `formula` is already a compiled design. For backward compatibility, the second
#'   positional argument after a compiled design may be the target term name.
#' @param distance A `dist` object or full distance matrix, used with a
#'   one-sided formula.
#' @param term Optional name of the sole categorical fixed main-effect term.
#'   It is inferred when the model contains exactly one fixed term.
#' @param type Group spatial `"median"` or `"centroid"`.
#' @param bias_adjust Apply the small-sample group-size adjustment.
#' @param permutations A non-negative integer or `"all"`.
#' @param seed Optional permutation seed.
#' @param keep_permutations Retain permutation statistics.
#' @param max_exact Maximum exact permutation space.
#' @param tol Numerical tolerance.
#' @param design Backward-compatible named alias for a compiled Design Object;
#'   use `formula = design` or pass the object as the first argument.
#' @return An object of class `permanova_dispersion`.
#' @examples
#' data(dune, dune.env, package = "vegan")
#' selected <- dune.env$Management %in% c("NM", "SF")
#' env <- droplevels(dune.env[selected, , drop = FALSE])
#' bray <- vegan::vegdist(dune[selected, ], method = "bray")
#' permanova_dispersion(
#'   bray ~ Management, env, permutations = 19, seed = 1
#' )
#' @export
permanova_dispersion <- function(
    formula = NULL, data = NULL, distance = NULL, term = NULL,
    type = c("median", "centroid"), bias_adjust = FALSE,
    permutations = 999L, seed = NULL, keep_permutations = FALSE,
    max_exact = 100000L, tol = 1e-10, design = NULL) {
  if (!is.null(design)) {
    if (!is.null(formula)) {
      stop("Supply the compiled design once, not through both `formula` and `design`.",
           call. = FALSE)
    }
    formula <- design
  }
  if (is.null(formula)) {
    stop("Supply a formula or a compiled `permanova_design` object.",
         call. = FALSE)
  }

  if (inherits(formula, "permanova_design")) {
    # Preserve the original `permanova_dispersion(design, "term")` call.
    if (is.null(term) && is.character(data) && length(data) == 1L &&
        !is.na(data) && nzchar(data)) {
      term <- data
      data <- NULL
    }
    if (!is.null(data)) {
      stop("`data` is not used when the first argument is a compiled design.",
           call. = FALSE)
    }
    if (!is.null(distance)) {
      stop("`distance` is already stored in the compiled design.", call. = FALSE)
    }
    compiled <- formula
  } else if (inherits(formula, "formula")) {
    model_formula <- formula
    if (length(formula) == 3L) {
      if (!is.null(distance)) {
        stop("Supply the dissimilarity once: on the formula's left-hand side or through `distance`.",
             call. = FALSE)
      }
      distance <- eval(
        formula[[2L]], envir = data, enclos = environment(formula)
      )
      model_formula <- stats::as.formula(
        call("~", formula[[3L]]), env = environment(formula)
      )
    }
    if (is.null(distance)) {
      stop(
        "Supply the dissimilarity: `permanova_dispersion(D ~ group, data)` or ",
        "`permanova_dispersion(~ group, data, distance = D)`.",
        call. = FALSE
      )
    }
    if (is.null(data)) {
      stop("`data` is required when `formula` is a model formula.", call. = FALSE)
    }
    compiled <- compile_permanova_design(
      model_formula, data, distance, tol = tol
    )
  } else {
    stop("`formula` must be a model formula or a compiled `permanova_design` object.",
         call. = FALSE)
  }

  if (is.null(term)) {
    if (length(compiled$fixed_effects) != 1L) {
      stop(
        "`term` can be inferred only when the model contains exactly one fixed term.",
        call. = FALSE
      )
    }
    term <- compiled$fixed_effects[[1L]]
  }

  type <- match.arg(type)
  target <- .validate_dispersion_target(compiled, term)
  geometry <- .dispersion_geometry(
    compiled$distance, target$group, type = type,
    bias_adjust = bias_adjust, tol = tol
  )
  response_design <- .dispersion_response_design(
    compiled, geometry$distances, tol = tol
  )
  if (!is.null(seed)) set.seed(seed)
  setup <- .dispersion_plan(
    compiled, response_design, permutations = permutations,
    max_exact = max_exact, tol = tol
  )
  if (nrow(setup$plan$orders) > 0L &&
      is.finite(setup$plan$attainable) && setup$plan$attainable < 100) {
    warning(
      "Only ", setup$plan$attainable,
      " dispersion-test permutations are attainable; interpret p-values cautiously.",
      call. = FALSE
    )
  }
  statistics <- .dispersion_statistics(
    geometry$distances, response_design, setup, tol = tol
  )
  notes <- geometry$warnings
  if (geometry$negative_count) {
    notes <- c(
      notes,
      "Negative principal-coordinate eigenvalues were retained in pseudo-Euclidean distances to group centers."
    )
  }
  exact <- setup$plan$exact && nrow(setup$plan$orders) > 0L
  specification <- list(
    term = term,
    hypothesis = paste0(
      "mean distance to the group ", type,
      " is equal across levels of `", term, "`"
    ),
    response = "distance to target-group center",
    center = type, bias_adjust = isTRUE(bias_adjust),
    numerator = list(
      ss = statistics$numerator_ss, df = setup$numerator$df,
      ms = statistics$numerator_ms, projector = setup$numerator$H
    ),
    denominator = list(
      term = setup$denominator$term, ss = statistics$denominator_ss,
      df = setup$denominator$df, ms = statistics$denominator_ms,
      projector = setup$denominator$H
    ),
    symbolic_ems = list(
      numerator = setup$numerator$ems,
      denominator = setup$denominator$ems,
      exact_match = max(abs(
        setup$numerator$ems - setup$denominator$ems
      )) <= sqrt(tol)
    ),
    statistic = statistics$observed, p_value = statistics$p_value,
    engine = setup$engine,
    exchangeable_units = setup$exchangeable_units,
    permutation_strata = setup$permutation_strata,
    permutation_algorithm = setup$algorithm,
    attainable_permutations = setup$plan$attainable,
    permutations_used = nrow(setup$plan$orders),
    unique_permutations_used = if (nrow(setup$plan$orders)) {
      nrow(unique(setup$plan$orders))
    } else 0L,
    exact_permutation = exact, permutation_seed = seed,
    approximation_status = if (!nrow(setup$plan$orders)) {
      "not tested"
    } else if (exact) {
      "exhaustive restricted full-model residual permutation"
    } else {
      "Monte Carlo restricted full-model residual permutation"
    },
    warnings = notes,
    permutation_statistics = if (keep_permutations) {
      statistics$statistics
    } else NULL
  )
  class(specification) <- "permanova_test_spec"

  table <- data.frame(
    term = c(term, setup$denominator$term),
    df = c(setup$numerator$df, setup$denominator$df),
    ss = c(statistics$numerator_ss, statistics$denominator_ss),
    ms = c(statistics$numerator_ms, statistics$denominator_ms),
    pseudo_f = c(statistics$observed, NA_real_),
    p_perm = c(statistics$p_value, NA_real_),
    row.names = NULL
  )
  out <- list(
    call = match.call(), design = compiled, term = term,
    engine = setup$engine, table = table, test = specification,
    distances = geometry$distances,
    group = geometry$group, group_distances = geometry$group_distances,
    centers = geometry$centers, coordinates = geometry$coordinates,
    eigenvalues = geometry$eigenvalues,
    retained_eigenvalues = geometry$retained_eigenvalues,
    geometry = geometry[c(
      "type", "bias_adjust", "positive_axes", "negative_axes",
      "negative_count", "zero_clamped", "warnings"
    )],
    permutation = list(
      algorithm = setup$algorithm, exact = exact,
      attainable = setup$plan$attainable,
      used = nrow(setup$plan$orders), seed = seed
    ),
    interpretation = paste(
      "This tests multivariate dispersion, not centroid location.",
      "A significant result can complicate location-test interpretation;",
      "a non-significant result does not prove exchangeability."
    )
  )
  class(out) <- "permanova_dispersion"
  out
}

#' @export
print.permanova_dispersion <- function(x, ...) {
  cat("Design-aware multivariate dispersion test\n")
  cat("Term:", x$term, " Center:", x$geometry$type, "\n")
  cat("Engine:", x$engine, "\n")
  cat("Permutation:", x$permutation$algorithm, "\n")
  print(x$table, row.names = FALSE)
  if (length(x$geometry$warnings)) {
    cat("Warnings:", paste(x$geometry$warnings, collapse = " "), "\n")
  }
  invisible(x)
}
