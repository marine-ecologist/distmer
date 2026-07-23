.hybrid_random_terms <- function(design, context) {
  if (!inherits(design, "permanova_design") || is.null(design$G)) {
    stop("`design` must be a compiled Design Object containing `distance`.",
         call. = FALSE)
  }
  fixed <- .one_categorical_fixed(design, context)
  if (length(design$random_terms) != 4L) {
    stop(
      paste0(
        context, " requires region, region/site, treatment:region, and ",
        "treatment:region:site random intercepts."
      ),
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
      paste0(
        context, " requires one region, one nested site term, and two ",
        "explicit treatment interactions."
      ),
      call. = FALSE
    )
  }
  region_term <- main_terms[[1L]]
  site_term <- nested_terms[[1L]]
  if (!identical(site_term$parent, region_term$label) ||
      length(region_term$variables) != 1L ||
      fixed$name %in% site_term$variables) {
    stop("The hybrid hierarchy must be one region/site containment pair distinct from treatment.",
         call. = FALSE)
  }
  region_name <- region_term$variables[[1L]]
  site_variables <- setdiff(site_term$variables, region_name)
  if (length(site_variables) != 1L) {
    stop("The hybrid hierarchy requires exactly one site variable within region.",
         call. = FALSE)
  }
  site_name <- site_variables[[1L]]
  interaction_variables <- lapply(interaction_terms, `[[`, "variables")
  region_index <- which(vapply(
    interaction_variables,
    function(value) setequal(value, c(fixed$name, region_name)), logical(1)
  ))
  site_index <- which(vapply(
    interaction_variables,
    function(value) setequal(value, c(fixed$name, region_name, site_name)),
    logical(1)
  ))
  if (length(region_index) != 1L || length(site_index) != 1L) {
    stop(
      "Declare explicit `(1 | treatment:region)` and `(1 | treatment:region:site)` covariance components.",
      call. = FALSE
    )
  }
  list(
    fixed = fixed,
    region_name = region_name, site_name = site_name,
    region_term = region_term$label, site_term = site_term$label,
    treatment_region_term = interaction_terms[[region_index]]$label,
    treatment_site_term = interaction_terms[[site_index]]$label
  )
}

.aggregate_hybrid_profiles <- function(
    design, treatment, region, site, treatment_site, treatment_region,
    tol = 1e-10) {
  cell_blocks <- unname(split(seq_len(design$n), treatment_site))
  cell_mapping <- do.call(rbind, lapply(seq_along(cell_blocks), function(index) {
    observations <- cell_blocks[[index]]
    values <- lapply(
      list(treatment, region, site),
      function(value) unique(as.character(value[observations]))
    )
    if (any(lengths(values) != 1L)) {
      stop("Every treatment-by-site cell must map to one treatment, region, and nested site.",
           call. = FALSE)
    }
    data.frame(
      cell = index, treatment = values[[1L]], region = values[[2L]],
      site = values[[3L]], size = length(observations),
      stringsAsFactors = FALSE
    )
  }))
  if (anyDuplicated(cell_mapping[c("treatment", "region", "site")])) {
    stop("Every observed treatment-by-region-by-site combination must form one cell.",
         call. = FALSE)
  }

  site_mapping <- unique(data.frame(
    site = as.character(site), region = as.character(region),
    stringsAsFactors = FALSE
  ))
  if (anyDuplicated(site_mapping$site)) {
    stop("Every nested site must belong to exactly one region.", call. = FALSE)
  }
  sites_per_region <- table(site_mapping$region)
  if (any(sites_per_region < 2L)) {
    stop("Unbalanced hybrid inference requires at least two sites in every region.",
         call. = FALSE)
  }

  profile_blocks <- unname(split(seq_len(design$n), treatment_region))
  profile_mapping <- do.call(rbind, lapply(
    seq_along(profile_blocks), function(index) {
      observations <- profile_blocks[[index]]
      treatment_value <- unique(as.character(treatment[observations]))
      region_value <- unique(as.character(region[observations]))
      if (length(treatment_value) != 1L || length(region_value) != 1L) {
        stop("Every treatment-by-region profile must map to one treatment and region.",
             call. = FALSE)
      }
      cells <- cell_mapping[
        cell_mapping$treatment == treatment_value &
          cell_mapping$region == region_value, , drop = FALSE
      ]
      expected_sites <- sort(site_mapping$site[
        site_mapping$region == region_value
      ])
      observed_sites <- sort(cells$site)
      if (!identical(observed_sites, expected_sites)) {
        stop(
          paste0(
            "Every observed treatment-by-region profile must contain the ",
            "complete site set for its region; partial profiles are unsupported."
          ),
          call. = FALSE
        )
      }
      cells <- cells[match(expected_sites, cells$site), , drop = FALSE]
      data.frame(
        profile = index, treatment = treatment_value,
        region = region_value, sites = nrow(cells),
        observations = sum(cells$size),
        topology = paste(cells$site, cells$size, sep = "=", collapse = ";"),
        stringsAsFactors = FALSE
      )
    }
  ))
  if (anyDuplicated(profile_mapping[c("treatment", "region")])) {
    stop("Every observed treatment-by-region combination must form one profile.",
         call. = FALSE)
  }

  membership <- matrix(
    0, nrow = nrow(profile_mapping), ncol = design$n,
    dimnames = list(NULL, seq_len(design$n))
  )
  for (profile_index in seq_len(nrow(profile_mapping))) {
    profile <- profile_mapping[profile_index, ]
    cells <- cell_mapping[
      cell_mapping$treatment == profile$treatment &
        cell_mapping$region == profile$region, , drop = FALSE
    ]
    for (cell_index in cells$cell) {
      observations <- cell_blocks[[cell_index]]
      membership[profile_index, observations] <-
        1 / (nrow(cells) * length(observations))
    }
  }
  if (max(abs(rowSums(membership) - 1)) > tol) {
    stop("Hybrid profile centroid weights do not sum to one.", call. = FALSE)
  }
  number_profiles <- nrow(membership)
  center <- diag(number_profiles) -
    matrix(1 / number_profiles, number_profiles, number_profiles)
  G <- center %*% (membership %*% design$G %*% t(membership)) %*% center
  G <- (G + t(G)) / 2
  G[abs(G) < tol] <- 0
  list(
    G = G, membership = membership, blocks = profile_blocks,
    sizes = profile_mapping$observations,
    representatives = vapply(profile_blocks, `[[`, integer(1), 1L),
    cell_blocks = cell_blocks, cell_mapping = cell_mapping,
    profile_mapping = profile_mapping, sites_per_region = sites_per_region
  )
}

.validate_unbalanced_hybrid <- function(design, tol = 1e-10) {
  terms <- .hybrid_random_terms(design, "Unbalanced hybrid inference")
  treatment <- droplevels(terms$fixed$value)
  region <- droplevels(design$groups[[terms$region_term]])
  site <- droplevels(design$groups[[terms$site_term]])
  treatment_region <- droplevels(
    design$groups[[terms$treatment_region_term]]
  )
  treatment_site <- droplevels(
    design$groups[[terms$treatment_site_term]]
  )
  if (nlevels(region) < 2L) {
    stop("Unbalanced hybrid inference requires at least two regions.",
         call. = FALSE)
  }
  aggregated <- .aggregate_hybrid_profiles(
    design, treatment, region, site, treatment_site, treatment_region, tol
  )
  mapping <- aggregated$profile_mapping
  profile_treatment <- factor(mapping$treatment, levels = levels(treatment))
  profile_region <- factor(mapping$region, levels = levels(region))
  incidence <- table(profile_region, profile_treatment)
  complete_incidence <- nrow(mapping) == nlevels(region) * nlevels(treatment) &&
    length(incidence) == nlevels(region) * nlevels(treatment) &&
    all(incidence == 1L)
  equal_sites <- length(unique(as.integer(aggregated$sites_per_region))) == 1L
  equal_cells <- length(unique(aggregated$cell_mapping$size)) == 1L
  if (complete_incidence && equal_sites && equal_cells) {
    stop("This hierarchy is balanced; use `permanova_hybrid()` for exact EMS inference.",
         call. = FALSE)
  }

  full_X <- stats::model.matrix(~ profile_region + profile_treatment)
  reduced_X <- stats::model.matrix(~ profile_region)
  full <- projection_matrix(full_X, tol = tol)
  reduced <- projection_matrix(reduced_X, tol = tol)
  numerator_df <- full$rank - reduced$rank
  denominator_df <- nrow(full_X) - full$rank
  if (full$rank != ncol(full_X) ||
      numerator_df != nlevels(profile_treatment) - 1L ||
      denominator_df <= 0L) {
    stop(
      paste0(
        "The treatment-by-region incidence is disconnected, aliases a ",
        "treatment contrast, or leaves no profile-centroid residual df."
      ),
      call. = FALSE
    )
  }
  exchangeability_strata <- interaction(
    profile_region, factor(mapping$topology), drop = TRUE, lex.order = TRUE
  )
  treatment_columns <- stats::model.matrix(~ 0 + profile_treatment)
  if (.stratified_treatment_rank(
        treatment_columns, exchangeability_strata, tol
      ) < numerator_df) {
    stop(
      paste0(
        "Region-by-profile-topology strata do not span every treatment ",
        "contrast; intact profile permutation cannot test the full hypothesis."
      ),
      call. = FALSE
    )
  }

  list(
    fixed = terms$fixed,
    region_name = terms$region_name, site_name = terms$site_name,
    region_term = terms$region_term, site_term = terms$site_term,
    treatment_region_term = terms$treatment_region_term,
    treatment_site_term = terms$treatment_site_term,
    treatment = profile_treatment, region = profile_region,
    aggregated = aggregated, incidence = incidence,
    exchangeability_strata = exchangeability_strata,
    H_reduced = reduced$H, H_numerator = full$H - reduced$H,
    numerator_df = numerator_df,
    H_denominator = diag(nrow(full_X)) - full$H,
    denominator_df = denominator_df,
    structure = list(
      type = "unbalanced crossed-nested treatment by region/site hybrid",
      fixed = terms$fixed, aggregated = aggregated,
      treatment = profile_treatment, model_strata = profile_region,
      permutation_strata = exchangeability_strata,
      strata_name = "region crossed with exact site-specific cell-count topology",
      unit_name = terms$treatment_region_term,
      estimand = paste(
        "equally weighted complete-site treatment-by-region profile centroids",
        "with additive region adjustment"
      ),
      denominator_name = "profile-centroid Residual after region and treatment",
      permutation_algorithm = paste(
        "contrast-specific Freedman-Lane permutation of intact profile",
        "centroids within region-by-topology strata"
      )
    )
  )
}

#' Unbalanced crossed-nested hybrid PERMANOVA
#'
#' Tests one treatment crossed with a region/site hierarchy when regions have
#' unequal site counts, treatment-by-region profiles are missing, or replicated
#' treatment-by-site cells have unequal sizes. Each observed profile must still
#' contain every site in its region. A profile is the equal-weight mean of its
#' site-cell centroids; observed profiles are adjusted for region and permuted
#' only within region and identical site-specific cell-count topology.
#'
#' @param design Canonical unbalanced hybrid Design Object with a distance.
#' @param permutations A non-negative integer or `"all"`.
#' @param seed Optional permutation seed.
#' @param keep_permutations Retain permutation statistics.
#' @param max_exact Maximum exact permutation space.
#' @param tol Numerical tolerance.
#' @return An object of class `permanova_unbalanced_hybrid`.
#' @export
permanova_unbalanced_hybrid <- function(
    design, permutations = 999L, seed = NULL,
    keep_permutations = FALSE, max_exact = 100000L, tol = 1e-10) {
  structure <- .validate_unbalanced_hybrid(design, tol = tol)
  if (!is.null(seed)) set.seed(seed)
  plan <- .stratified_permutation_plan(
    structure$exchangeability_strata, permutations, max_exact
  )
  result <- .unit_raw_permutation(
    structure$aggregated$G, structure$H_numerator,
    structure$numerator_df, structure$H_denominator,
    structure$denominator_df, plan, tol
  )
  notes <- paste(
    "Inference is empirical on profile centroids; no unbalanced parametric",
    "EMS or variance-component claim is made."
  )
  if (!isTRUE(design$spectral$euclidean)) {
    notes <- c(
      notes,
      "Empirical inference retained the complete non-Euclidean Gower spectrum."
    )
  }
  test <- list(
    term = structure$fixed$name,
    hypothesis = paste0(
      "complete-site profile centroids are equal across `",
      structure$fixed$name, "` levels after region adjustment"
    ),
    estimand = structure$structure$estimand,
    numerator = list(
      ss = result$numerator_ss, df = structure$numerator_df,
      ms = result$numerator_ms, projector = structure$H_numerator
    ),
    denominator = list(
      expression = "MS[profile-centroid Residual after region and treatment]",
      ss = result$denominator_ss, df = structure$denominator_df,
      ms = result$denominator_ms, projector = structure$H_denominator
    ),
    statistic = result$observed, p_value = result$p_value,
    exchangeable_units = structure$treatment_region_term,
    descendant_units = c(
      structure$treatment_site_term, structure$site_term, "observations"
    ),
    permutation_strata = structure$structure$strata_name,
    permutation_algorithm = paste(
      "raw permutation of intact equal-site-weight profile centroids within",
      "region-by-exact-topology strata"
    ),
    attainable_permutations = plan$attainable,
    permutations_used = nrow(plan$orders),
    unique_permutations_used = if (nrow(plan$orders)) {
      nrow(unique(plan$orders))
    } else 0L,
    exact_permutation = plan$exact && nrow(plan$orders) > 0L,
    approximation_status = if (!nrow(plan$orders)) {
      "not tested"
    } else if (plan$exact) {
      "exhaustive region-by-topology profile permutation"
    } else "Monte Carlo region-by-topology profile permutation",
    warnings = notes,
    permutation_statistics = if (keep_permutations) {
      result$statistics
    } else NULL
  )
  class(test) <- "permanova_test_spec"
  table <- data.frame(
    term = c(structure$fixed$name, "Profile-centroid Residual"),
    df = c(structure$numerator_df, structure$denominator_df),
    ss = c(result$numerator_ss, result$denominator_ss),
    ms = c(result$numerator_ms, result$denominator_ms),
    pseudo_f = c(result$observed, NA_real_),
    p_perm = c(result$p_value, NA_real_)
  )
  out <- list(
    call = match.call(), design = design, table = table, test = test,
    structure = list(
      treatment = structure$fixed$name, region = structure$region_name,
      site = structure$site_name,
      observed_profiles = nrow(structure$aggregated$profile_mapping),
      sites_per_region = structure$aggregated$sites_per_region,
      profile_topologies = structure$aggregated$profile_mapping$topology,
      incidence = structure$incidence
    ),
    estimand = test$estimand,
    supported_design = paste(
      "one treatment crossed with one unbalanced region/site hierarchy;",
      "complete observed site profiles and empirical topology-stratified inference"
    ),
    permutation = list(
      algorithm = test$permutation_algorithm,
      exact = test$exact_permutation, attainable = plan$attainable,
      used = nrow(plan$orders), seed = seed
    )
  )
  class(out) <- "permanova_unbalanced_hybrid"
  out
}

#' @export
print.permanova_unbalanced_hybrid <- function(x, ...) {
  cat("Unbalanced crossed-nested hybrid PERMANOVA\n")
  cat("Estimand:", x$estimand, "\n")
  cat("Permutation:", x$permutation$algorithm, "\n")
  print(x$table, row.names = FALSE)
  invisible(x)
}
