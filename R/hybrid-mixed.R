.validate_hybrid_design <- function(design, tol = 1e-10) {
  if (!inherits(design, "permanova_design")) {
    stop("`design` must be produced by `compile_permanova_design()`.",
         call. = FALSE)
  }
  if (is.null(design$G)) {
    stop("The Design Object must include `distance`.", call. = FALSE)
  }
  fixed_variables <- all.vars(design$fixed_formula)
  if (length(design$fixed_effects) != 1L || length(fixed_variables) != 1L ||
      !identical(.canonical_term(design$fixed_effects[[1L]]),
                 fixed_variables[[1L]])) {
    stop("Hybrid inference supports exactly one categorical fixed treatment factor.",
         call. = FALSE)
  }
  treatment_name <- fixed_variables[[1L]]
  if (!is.factor(design$data[[treatment_name]])) {
    stop("Hybrid inference requires a categorical fixed treatment with at least two levels.",
         call. = FALSE)
  }
  treatment <- droplevels(design$data[[treatment_name]])
  if (nlevels(treatment) < 2L) {
    stop("Hybrid inference requires a categorical fixed treatment with at least two levels.",
         call. = FALSE)
  }
  if (length(design$estimability$aliased) ||
      any(!design$estimability$terms)) {
    stop("The hybrid fixed-effect hypothesis is rank deficient or aliased.",
         call. = FALSE)
  }
  if (length(design$random_terms) != 4L) {
    stop(
      "Hybrid inference requires region, region/site, treatment:region, and treatment:region:site random intercepts.",
      call. = FALSE
    )
  }

  origins <- vapply(design$random_terms, `[[`, character(1), "origin")
  main_terms <- design$random_terms[origins == "main"]
  nested_terms <- design$random_terms[origins == "nesting"]
  interaction_terms <- design$random_terms[origins == "interaction"]
  if (length(main_terms) != 1L || length(nested_terms) != 1L ||
      length(interaction_terms) != 2L) {
    stop(
      "Hybrid random provenance must contain one region, one nested site term, and two explicit treatment interactions.",
      call. = FALSE
    )
  }
  region_term <- main_terms[[1L]]
  nested_term <- nested_terms[[1L]]
  if (!identical(nested_term$parent, region_term$label) ||
      length(region_term$variables) != 1L ||
      treatment_name %in% nested_term$variables) {
    stop("The hybrid hierarchy must be one region/site containment pair distinct from treatment.",
         call. = FALSE)
  }
  region_name <- region_term$variables[[1L]]
  child_variables <- setdiff(nested_term$variables, region_name)
  if (length(child_variables) != 1L) {
    stop("The Phase 9 hybrid hierarchy requires exactly one site variable within region.",
         call. = FALSE)
  }
  site_name <- child_variables[[1L]]
  interaction_variables <- lapply(interaction_terms, `[[`, "variables")
  region_interaction_index <- which(vapply(
    interaction_variables,
    function(x) setequal(x, c(treatment_name, region_name)), logical(1)
  ))
  site_interaction_index <- which(vapply(
    interaction_variables,
    function(x) setequal(x, c(treatment_name, region_name, site_name)),
    logical(1)
  ))
  if (length(region_interaction_index) != 1L ||
      length(site_interaction_index) != 1L) {
    stop(
      "Declare explicit `(1 | treatment:region)` and `(1 | treatment:region:site)` covariance components.",
      call. = FALSE
    )
  }
  treatment_region_term <- interaction_terms[[region_interaction_index]]
  treatment_site_term <- interaction_terms[[site_interaction_index]]

  region <- droplevels(design$groups[[region_term$label]])
  site_unit <- droplevels(design$groups[[nested_term$label]])
  treatment_region <- droplevels(
    design$groups[[treatment_region_term$label]]
  )
  treatment_site <- droplevels(
    design$groups[[treatment_site_term$label]]
  )
  if (nlevels(region) < 2L) {
    stop("Hybrid inference requires at least two regions.", call. = FALSE)
  }

  site_mapping <- unique(data.frame(
    site = as.integer(site_unit), region = as.integer(region)
  ))
  if (any(table(site_mapping$site) != 1L)) {
    stop("Every nested site must belong to exactly one region.", call. = FALSE)
  }
  sites_per_region <- table(site_mapping$region)
  if (length(unique(as.integer(sites_per_region))) != 1L ||
      sites_per_region[[1L]] < 2L) {
    stop("Balanced hybrid inference requires at least two sites in every region.",
         call. = FALSE)
  }

  cell_mapping <- unique(data.frame(
    cell = as.integer(treatment_site),
    treatment = as.integer(treatment),
    region = as.integer(region), site = as.integer(site_unit)
  ))
  if (any(table(cell_mapping$cell) != 1L)) {
    stop("Every treatment-by-site cell must map to one treatment and nested site.",
         call. = FALSE)
  }
  incidence <- table(cell_mapping$site, cell_mapping$treatment)
  if (any(incidence != 1L)) {
    stop("Every treatment must occur exactly once within every nested site.",
         call. = FALSE)
  }
  cell_blocks <- unname(split(seq_len(design$n), treatment_site))
  if (length(unique(lengths(cell_blocks))) != 1L ||
      lengths(cell_blocks)[[1L]] < 2L) {
    stop(
      "Balanced hybrid inference requires equal replication of at least two observations per treatment-by-site cell.",
      call. = FALSE
    )
  }

  profile_mapping <- unique(data.frame(
    profile = as.integer(treatment_region),
    treatment = as.integer(treatment), region = as.integer(region)
  ))
  if (any(table(profile_mapping$profile) != 1L)) {
    stop("Every treatment profile must belong to one treatment and one region.",
         call. = FALSE)
  }
  profile_incidence <- table(
    profile_mapping$region, profile_mapping$treatment
  )
  if (any(profile_incidence != 1L)) {
    stop("Every treatment must have one complete profile within every region.",
         call. = FALSE)
  }
  profile_blocks <- unname(split(seq_len(design$n), treatment_region))
  if (length(unique(lengths(profile_blocks))) != 1L) {
    stop("Every treatment-by-region profile must contain equal descendants.",
         call. = FALSE)
  }
  profile_mapping <- profile_mapping[
    order(profile_mapping$region, profile_mapping$treatment), , drop = FALSE
  ]
  profiles_by_region <- unname(split(
    profile_mapping$profile, profile_mapping$region
  ))

  P0 <- projection_matrix(matrix(1, design$n, 1L), tol = tol)
  PT <- design$projections$fixed
  PR <- design$projections$random[[region_term$label]]
  PRS <- design$projections$random[[nested_term$label]]
  PTR <- design$projections$random[[treatment_region_term$label]]
  PTRS <- design$projections$random[[treatment_site_term$label]]
  projectors <- list(
    treatment = PT$H - P0$H,
    region = PR$H - P0$H,
    `site within region` = PRS$H - PR$H,
    `treatment:region` = PTR$H - PT$H - PR$H + P0$H,
    `treatment:site within region` =
      PTRS$H - PTR$H - PRS$H + PR$H,
    Residual = diag(design$n) - PTRS$H
  )
  dfs <- c(
    treatment = PT$rank - P0$rank,
    region = PR$rank - P0$rank,
    `site within region` = PRS$rank - PR$rank,
    `treatment:region` = PTR$rank - PT$rank - PR$rank + P0$rank,
    `treatment:site within region` =
      PTRS$rank - PTR$rank - PRS$rank + PR$rank,
    Residual = design$n - PTRS$rank
  )
  if (any(dfs <= 0L)) {
    stop("Every hybrid ANOVA stratum must retain positive degrees of freedom.",
         call. = FALSE)
  }
  projector_tolerance <- max(100 * tol, sqrt(.Machine$double.eps))
  idempotence <- vapply(projectors, function(H) {
    max(abs(H %*% H - H))
  }, numeric(1))
  pair_indices <- utils::combn(seq_along(projectors), 2L, simplify = FALSE)
  orthogonality <- max(vapply(pair_indices, function(pair) {
    max(abs(projectors[[pair[[1L]]]] %*% projectors[[pair[[2L]]]]))
  }, numeric(1)))
  if (max(idempotence) > projector_tolerance ||
      orthogonality > projector_tolerance) {
    stop(
      "Hybrid factor projections are not a balanced orthogonal crossed-nested decomposition.",
      call. = FALSE
    )
  }

  components <- c(
    "residual", region_term$label, nested_term$label,
    treatment_region_term$label, treatment_site_term$label
  )
  covariance_bases <- c(
    list(residual = diag(design$n)),
    design$Z[c(
      region_term$label, nested_term$label,
      treatment_region_term$label, treatment_site_term$label
    )]
  )
  covariance_bases[-1L] <- lapply(covariance_bases[-1L], tcrossprod)
  covariance_bases <- covariance_bases[components]
  strata <- setNames(lapply(seq_along(projectors), function(index) {
    H <- projectors[[index]]
    df <- unname(dfs[[index]])
    ss <- trace_hgh(H, design$G)
    list(
      term = names(projectors)[[index]], H = H, df = df,
      ss = ss, ms = ss / df,
      ems = setNames(vapply(covariance_bases, function(V) {
        trace_product(H, V) / df
      }, numeric(1)), components)
    )
  }), names(projectors))
  if (max(abs(
    strata$treatment$ems - strata$`treatment:region`$ems
  )) > sqrt(tol)) {
    stop("Treatment and treatment-by-region strata do not have matching null EMS.",
         call. = FALSE)
  }

  list(
    treatment_name = treatment_name, region_name = region_name,
    site_name = site_name, treatment = treatment, region = region,
    site_unit = site_unit, treatment_region = treatment_region,
    treatment_site = treatment_site,
    region_term = region_term$label, site_term = nested_term$label,
    treatment_region_term = treatment_region_term$label,
    treatment_site_term = treatment_site_term$label,
    treatments = nlevels(treatment), regions = nlevels(region),
    sites_per_region = unname(sites_per_region[[1L]]),
    replicates = lengths(cell_blocks)[[1L]],
    cell_blocks = cell_blocks, profile_blocks = profile_blocks,
    profiles_by_region = profiles_by_region,
    profile_mapping = profile_mapping,
    components = components, covariance_bases = covariance_bases,
    strata = strata, idempotence_error = max(idempotence),
    orthogonality_error = orthogonality
  )
}

#' Derive EMS for a balanced crossed-nested hybrid
#'
#' @param design Canonical Phase 9 hybrid Design Object.
#' @param tol Numerical tolerance.
#' @return A machine-readable treatment and denominator EMS specification.
#' @export
derive_hybrid_ems <- function(design, tol = 1e-10) {
  structure <- .validate_hybrid_design(design, tol = tol)
  numerator <- structure$strata$treatment
  denominator <- structure$strata$`treatment:region`
  list(
    term = structure$treatment_name,
    hypothesis = paste0(
      "all estimable coefficients for `", structure$treatment_name,
      "` are zero"
    ),
    numerator_df = numerator$df,
    numerator_under_null = numerator$ems,
    numerator_under_alternative = c(numerator$ems, fixed_effect = 1),
    denominator = denominator$ems,
    denominator_term = denominator$term,
    denominator_df = denominator$df,
    exact_denominator = TRUE,
    derivation = paste0(
      "Treatment profiles cross regions and contain every nested site. ",
      "Trace-derived EMS for MS[", structure$treatment_name,
      "] and MS[", denominator$term,
      "] match for residual, hierarchy, and declared treatment-by-hierarchy covariance components under H0."
    )
  )
}

.hybrid_permutation_plan <- function(structure, permutations, max_exact) {
  .crossed_block_permutation_plan(
    list(
      treatments = structure$treatments, blocks = structure$regions,
      cell_blocks = structure$profile_blocks,
      cells_by_block = structure$profiles_by_region,
      interaction_name = structure$treatment_region_term,
      block_name = structure$region_name
    ),
    permutations = permutations, max_exact = max_exact
  )
}

.hybrid_permutation_statistics <- function(design, structure, plan, tol) {
  numerator <- structure$strata$treatment
  denominator <- structure$strata$`treatment:region`
  if (abs(denominator$ms) <= tol) {
    stop("The treatment-by-region denominator mean square is numerically zero.",
         call. = FALSE)
  }
  observed <- numerator$ms / denominator$ms
  statistics <- numeric(nrow(plan$orders))
  for (row in seq_len(nrow(plan$orders))) {
    index <- .block_permutation_index(
      structure$profile_blocks, plan$orders[row, ]
    )
    G_permuted <- design$G[index, index, drop = FALSE]
    statistics[[row]] <-
      (trace_hgh(numerator$H, G_permuted) / numerator$df) /
      (trace_hgh(denominator$H, G_permuted) / denominator$df)
  }
  p_value <- NA_real_
  if (length(statistics)) {
    extreme <- sum(statistics >= observed - tol, na.rm = TRUE)
    p_value <- if (plan$exact) {
      extreme / length(statistics)
    } else {
      (extreme + 1) / (length(statistics) + 1)
    }
  }
  list(observed = observed, statistics = statistics, p_value = p_value)
}

#' Estimate balanced crossed-nested hybrid variance components
#'
#' @param design Canonical Phase 9 hybrid Design Object.
#' @param truncate Also report a separate zero-truncated sensitivity estimate.
#' @param tol Numerical tolerance.
#' @return A trace-derived hybrid variance-component specification.
#' @export
estimate_hybrid_variance_components <- function(
    design, truncate = FALSE, tol = 1e-10) {
  structure <- .validate_hybrid_design(design, tol = tol)
  selected_names <- c(
    "Residual", "treatment:site within region", "site within region",
    "treatment:region", "region"
  )
  selected <- structure$strata[selected_names]
  coefficient_matrix <- do.call(rbind, lapply(selected, `[[`, "ems"))
  observed_ms <- vapply(selected, `[[`, numeric(1), "ms")
  if (qr(coefficient_matrix, tol = tol)$rank < ncol(coefficient_matrix)) {
    stop("Hybrid variance-component EMS system is singular.", call. = FALSE)
  }
  estimates <- setNames(
    as.vector(solve(coefficient_matrix, observed_ms)), structure$components
  )
  list(
    unrestricted = estimates,
    truncated = if (truncate) pmax(estimates, 0) else NULL,
    coefficient_matrix = coefficient_matrix,
    observed_mean_squares = observed_ms,
    strata = selected_names,
    warnings = if (any(estimates < -tol)) {
      "Negative unrestricted hybrid variance components were retained."
    } else character(),
    method = "balanced trace-derived crossed-nested EMS equations"
  )
}

#' Balanced crossed-nested hybrid PERMANOVA
#'
#' Tests one fixed treatment crossed with regions and with sites nested within
#' regions. Required covariance components are region, region:site,
#' treatment:region, and treatment:region:site. Complete treatment-by-region
#' profiles permute within regions and carry every nested descendant intact.
#'
#' @param design Canonical Phase 9 hybrid Design Object.
#' @param permutations A non-negative integer or `"all"`.
#' @param seed Optional permutation seed.
#' @param keep_permutations Retain permutation statistics.
#' @param truncate Also report zero-truncated variance components.
#' @param max_exact Maximum exact permutation space.
#' @param tol Numerical tolerance.
#' @return An object of class `permanova_hybrid`.
#' @export
permanova_hybrid <- function(
    design, permutations = 999L, seed = NULL,
    keep_permutations = FALSE, truncate = FALSE,
    max_exact = 100000L, tol = 1e-10) {
  structure <- .validate_hybrid_design(design, tol = tol)
  ems <- derive_hybrid_ems(design, tol = tol)
  if (!is.null(seed)) set.seed(seed)
  plan <- .hybrid_permutation_plan(
    structure, permutations = permutations, max_exact = max_exact
  )
  if (nrow(plan$orders) > 0L && is.finite(plan$attainable) &&
      plan$attainable < 100) {
    warning(
      "Only ", plan$attainable,
      " treatment-profile permutations are attainable; interpret p-values cautiously.",
      call. = FALSE
    )
  }
  statistics <- .hybrid_permutation_statistics(
    design, structure, plan, tol = tol
  )
  numerator <- structure$strata$treatment
  denominator <- structure$strata$`treatment:region`
  notes <- character()
  if (!isTRUE(design$spectral$euclidean)) {
    notes <- c(
      notes,
      "Empirical inference retained the complete non-Euclidean Gower spectrum."
    )
  }
  if (numerator$ss < -tol || denominator$ss < -tol) {
    notes <- c(notes, "A negative distance-based sum of squares was retained.")
  }
  exact <- plan$exact && nrow(plan$orders) > 0L
  test <- list(
    term = structure$treatment_name, hypothesis = ems$hypothesis,
    numerator = numerator,
    denominator = list(
      expression = paste0("MS[", denominator$term, "]"),
      ss = denominator$ss, df = denominator$df, ms = denominator$ms,
      projector = denominator$H
    ),
    symbolic_ems = ems, statistic = statistics$observed,
    p_value = statistics$p_value,
    exchangeable_units = structure$treatment_region_term,
    permutation_strata = structure$region_name,
    descendant_units = c(
      structure$treatment_site_term, structure$site_term, "observations"
    ),
    permutation_algorithm = paste0(
      "raw permutation of complete `", structure$treatment_region_term,
      "` profiles within `", structure$region_name,
      "`, moving all nested sites and replicates intact"
    ),
    attainable_permutations = plan$attainable,
    permutations_used = nrow(plan$orders),
    unique_permutations_used = if (nrow(plan$orders)) {
      nrow(unique(plan$orders))
    } else 0L,
    exact_permutation = exact, permutation_seed = seed,
    approximation_status = if (!nrow(plan$orders)) {
      "not tested"
    } else if (exact) {
      "exhaustive within-region treatment-profile permutation"
    } else {
      "Monte Carlo within-region treatment-profile permutation"
    },
    warnings = notes,
    permutation_statistics = if (keep_permutations) {
      statistics$statistics
    } else NULL
  )
  class(test) <- "permanova_test_spec"
  table <- do.call(rbind, lapply(structure$strata, function(stratum) {
    data.frame(
      term = stratum$term, df = stratum$df, ss = stratum$ss,
      ms = stratum$ms,
      pseudo_f = if (stratum$term == "treatment") statistics$observed else NA_real_,
      p_perm = if (stratum$term == "treatment") statistics$p_value else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
  rownames(table) <- NULL
  out <- list(
    call = match.call(), design = design, table = table, test = test,
    ems = ems,
    variance_components = estimate_hybrid_variance_components(
      design, truncate = truncate, tol = tol
    ),
    structure = list(
      treatment = structure$treatment_name,
      region = structure$region_name, site = structure$site_name,
      treatments = structure$treatments, regions = structure$regions,
      sites_per_region = structure$sites_per_region,
      replicates = structure$replicates
    ),
    supported_design = paste(
      "balanced treatment crossed with one region/site hierarchy;",
      "declared treatment:region and treatment:region:site covariance components"
    ),
    permutation = list(
      algorithm = test$permutation_algorithm, exact = exact,
      attainable = plan$attainable, used = nrow(plan$orders), seed = seed
    )
  )
  class(out) <- "permanova_hybrid"
  out
}

#' @export
print.permanova_hybrid <- function(x, ...) {
  cat("Balanced crossed-nested hybrid PERMANOVA\n")
  cat("Supported design:", x$supported_design, "\n")
  cat("Permutation:", x$permutation$algorithm, "\n")
  print(x$table, row.names = FALSE)
  cat("Unrestricted variance components:\n")
  print(x$variance_components$unrestricted)
  invisible(x)
}
