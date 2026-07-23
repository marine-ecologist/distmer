.positive_integer <- function(x, name) {
  if (length(x) != 1L || !is.finite(x) || x < 1 || x != as.integer(x)) {
    stop("`", name, "` must be one positive integer.", call. = FALSE)
  }
  as.integer(x)
}

#' Simulate a balanced random-intercept community design
#'
#' @param groups_per_treatment Whole plots per treatment level.
#' @param replicates Observations per whole plot.
#' @param response_variables Number of response dimensions.
#' @param effect Euclidean location-shift magnitude.
#' @param cluster_sd Random-intercept standard deviation.
#' @param residual_sd Residual standard deviation in the control group.
#' @param dispersion_ratio Residual SD multiplier in the impact group.
#' @param seed Optional random seed.
#' @return A list with model data, response matrix, and Euclidean distance.
#' @export
simulate_balanced_mixed <- function(groups_per_treatment = 4L, replicates = 3L,
                                    response_variables = 1L, effect = 0,
                                    cluster_sd = 1, residual_sd = 1,
                                    dispersion_ratio = 1, seed = NULL) {
  groups_per_treatment <- .positive_integer(groups_per_treatment, "groups_per_treatment")
  replicates <- .positive_integer(replicates, "replicates")
  response_variables <- .positive_integer(response_variables, "response_variables")
  if (any(!is.finite(c(effect, cluster_sd, residual_sd, dispersion_ratio))) ||
      cluster_sd < 0 || residual_sd < 0 || dispersion_ratio < 0) {
    stop("Simulation effects and standard deviations must be finite and non-negative.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)

  groups <- 2L * groups_per_treatment
  data <- expand.grid(replicate = seq_len(replicates), site = factor(seq_len(groups)))
  data$treatment <- factor(
    ifelse(as.integer(data$site) <= groups_per_treatment, "control", "impact"),
    levels = c("control", "impact")
  )
  cluster_effect <- matrix(rnorm(groups * response_variables, sd = cluster_sd),
                           nrow = groups)
  residual_scale <- residual_sd * ifelse(data$treatment == "impact", dispersion_ratio, 1)
  residual <- matrix(rnorm(nrow(data) * response_variables), ncol = response_variables) *
    residual_scale
  direction <- rep(effect / sqrt(response_variables), response_variables)
  response <- cluster_effect[as.integer(data$site), , drop = FALSE] + residual +
    (data$treatment == "impact") * matrix(direction, nrow(data), response_variables, byrow = TRUE)
  colnames(response) <- paste0("community_", seq_len(response_variables))

  list(data = data, response = response, distance = stats::dist(response),
       parameters = list(
         groups_per_treatment = groups_per_treatment, replicates = replicates,
         response_variables = response_variables, effect = effect,
         cluster_sd = cluster_sd, residual_sd = residual_sd,
         dispersion_ratio = dispersion_ratio, seed = seed
       ))
}

#' Run mixed-model operating-characteristic simulations
#'
#' Evaluates null rejection, location-shift power, dispersion sensitivity,
#' denominator agreement with scalar ANOVA, failures, runtime, and numerical
#' stability. Rejection-rate uncertainty is reported as a Monte Carlo standard
#' error; small runs are smoke tests, not inferential validation.
#'
#' @param simulations Number of simulations per scenario.
#' @param permutations Permutations per fitted model.
#' @param alpha Rejection threshold.
#' @param groups_per_treatment,replicates Balanced design sizes.
#' @param effect Location shift for the power scenario.
#' @param dispersion_ratio Impact/control residual SD ratio for the dispersion scenario.
#' @param seed Simulation seed.
#' @param keep_results Retain simulation-level results.
#' @return A simulation-validation report.
#' @export
simulate_mixed_validation <- function(simulations = 100L, permutations = 199L,
                                      alpha = 0.05, groups_per_treatment = 4L,
                                      replicates = 3L, effect = 1,
                                      dispersion_ratio = 2, seed = 1L,
                                      keep_results = FALSE) {
  simulations <- .positive_integer(simulations, "simulations")
  groups_per_treatment <- .positive_integer(groups_per_treatment, "groups_per_treatment")
  replicates <- .positive_integer(replicates, "replicates")
  if (length(permutations) != 1L || permutations < 0 ||
      permutations != as.integer(permutations)) {
    stop("`permutations` must be one non-negative integer.", call. = FALSE)
  }
  if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must lie strictly between zero and one.", call. = FALSE)
  }

  scenarios <- data.frame(
    scenario = c("null", "location_shift", "dispersion_heterogeneity"),
    estimand = c("type_I_error", "power", "dispersion_robustness"),
    effect = c(0, effect, 0), dispersion = c(1, 1, dispersion_ratio),
    stringsAsFactors = FALSE
  )
  raw <- vector("list", nrow(scenarios))
  elapsed <- numeric(nrow(scenarios))

  for (scenario_index in seq_len(nrow(scenarios))) {
    start <- proc.time()[["elapsed"]]
    rows <- vector("list", simulations)
    for (iteration in seq_len(simulations)) {
      iteration_seed <- as.integer(seed + scenario_index * 100000L + iteration)
      rows[[iteration]] <- tryCatch({
        generated <- simulate_balanced_mixed(
          groups_per_treatment = groups_per_treatment, replicates = replicates,
          effect = scenarios$effect[[scenario_index]],
          dispersion_ratio = scenarios$dispersion[[scenario_index]], seed = iteration_seed
        )
        design <- compile_permanova_design(
          ~ treatment + (1 | site), generated$data, generated$distance
        )
        fit <- suppressWarnings(permanova_mixed(
          design, permutations = permutations, seed = iteration_seed
        ))
        site_mean <- aggregate(
          generated$response[, 1L],
          list(site = generated$data$site, treatment = generated$data$treatment), mean
        )
        names(site_mean)[[3L]] <- "response"
        reference_f <- unname(summary(aov(response ~ treatment, site_mean))[[1]][["F value"]][1])
        data.frame(
          iteration = iteration, p_value = fit$tests$treatment$p_value,
          pseudo_f = fit$tests$treatment$statistic,
          reference_f = reference_f, f_difference = fit$tests$treatment$statistic - reference_f,
          failed = FALSE, message = "", stringsAsFactors = FALSE
        )
      }, error = function(e) data.frame(
        iteration = iteration, p_value = NA_real_, pseudo_f = NA_real_,
        reference_f = NA_real_, f_difference = NA_real_, failed = TRUE,
        message = conditionMessage(e), stringsAsFactors = FALSE
      ))
    }
    elapsed[[scenario_index]] <- proc.time()[["elapsed"]] - start
    raw[[scenario_index]] <- do.call(rbind, rows)
    raw[[scenario_index]]$scenario <- scenarios$scenario[[scenario_index]]
  }

  summary <- do.call(rbind, lapply(seq_len(nrow(scenarios)), function(k) {
    values <- raw[[k]]
    valid <- !values$failed & is.finite(values$p_value)
    rate <- if (any(valid)) mean(values$p_value[valid] <= alpha) else NA_real_
    data.frame(
      scenario = scenarios$scenario[[k]], estimand = scenarios$estimand[[k]],
      simulations = simulations, successful = sum(valid), failures = sum(values$failed),
      rejection_rate = rate,
      monte_carlo_se = if (is.finite(rate)) sqrt(rate * (1 - rate) / sum(valid)) else NA_real_,
      mean_p_value = if (any(valid)) mean(values$p_value[valid]) else NA_real_,
      max_abs_denominator_f_difference = if (any(valid)) {
        max(abs(values$f_difference[valid]))
      } else NA_real_,
      elapsed_seconds = elapsed[[k]], stringsAsFactors = FALSE
    )
  }))
  rownames(summary) <- NULL
  out <- list(
    summary = summary, results = if (keep_results) do.call(rbind, raw) else NULL,
    settings = list(alpha = alpha, permutations = as.integer(permutations), seed = seed),
    interpretation = "Monte Carlo uncertainty must be considered; increase simulations for scientific validation."
  )
  class(out) <- "permanova_simulation_validation"
  out
}

#' @export
print.permanova_simulation_validation <- function(x, ...) {
  cat("PERMANOVA mixed-model simulation validation\n")
  print(x$summary, row.names = FALSE)
  cat(x$interpretation, "\n")
  invisible(x)
}

.binomial_summary <- function(rejected, simulations) {
  interval <- stats::binom.test(rejected, simulations)$conf.int
  c(
    rejection_rate = rejected / simulations,
    monte_carlo_se = sqrt((rejected / simulations) * (1 - rejected / simulations) / simulations),
    confidence_low = interval[[1L]], confidence_high = interval[[2L]]
  )
}

#' Validate fixed-by-random interaction Type I error
#'
#' Simulates a null treatment effect with genuine treatment-by-block random
#' variation. It compares the declared-interaction engine, the identical
#' classical interaction-denominator test, and an explicitly opted-in additive
#' model that omits the interaction variance.
#'
#' @param simulations Number of simulation replicates.
#' @param blocks Number of random blocks.
#' @param replicates Observations per treatment-by-block cell.
#' @param interaction_sd,block_sd,residual_sd Generating standard deviations.
#' @param alpha Rejection threshold.
#' @param seed Reproducibility seed.
#' @param keep_results Retain iteration-level results.
#' @return An interaction-validation report with binomial confidence intervals.
#' @export
simulate_crossed_interaction_validation <- function(
    simulations = 500L, blocks = 8L, replicates = 3L,
    interaction_sd = 2, block_sd = 1, residual_sd = 1,
    alpha = 0.05, seed = 1L, keep_results = FALSE) {
  simulations <- .positive_integer(simulations, "simulations")
  blocks <- .positive_integer(blocks, "blocks")
  replicates <- .positive_integer(replicates, "replicates")
  if (blocks < 3L || replicates < 2L) {
    stop("At least three blocks and two cell replicates are required.", call. = FALSE)
  }
  scales <- c(interaction_sd, block_sd, residual_sd)
  if (any(!is.finite(scales)) || any(scales < 0)) {
    stop("Generating standard deviations must be finite and non-negative.", call. = FALSE)
  }
  if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must lie strictly between zero and one.", call. = FALSE)
  }

  data <- expand.grid(
    replicate = factor(seq_len(replicates)),
    treatment = factor(c("control", "impact")),
    block = factor(seq_len(blocks))
  )
  rows <- vector("list", simulations)
  start <- proc.time()[["elapsed"]]
  for (iteration in seq_len(simulations)) {
    set.seed(as.integer(seed + iteration))
    block_effect <- rnorm(blocks, sd = block_sd)
    interaction_effect <- rnorm(2L * blocks, sd = interaction_sd)
    cell <- interaction(data$treatment, data$block, drop = TRUE)
    response <- block_effect[data$block] + interaction_effect[cell] +
      rnorm(nrow(data), sd = residual_sd)
    distance <- dist(response)

    declared_design <- compile_permanova_design(
      ~ treatment + (1 | block) + (1 | treatment:block), data, distance
    )
    declared <- permanova_crossed(declared_design, p_method = "satterthwaite")
    additive_design <- compile_permanova_design(
      ~ treatment + (1 | block), data, distance
    )
    additive <- permanova_crossed(
      additive_design, p_method = "satterthwaite", assume_additive = TRUE
    )

    reference <- summary(aov(response ~ treatment * block, data))[[1]]
    reference_rows <- trimws(rownames(reference))
    treatment_row <- match("treatment", reference_rows)
    interaction_row <- match("treatment:block", reference_rows)
    classical_f <- reference[["Mean Sq"]][[treatment_row]] /
      reference[["Mean Sq"]][[interaction_row]]
    classical_p <- stats::pf(
      classical_f, reference[["Df"]][[treatment_row]],
      reference[["Df"]][[interaction_row]], lower.tail = FALSE
    )
    rows[[iteration]] <- data.frame(
      iteration = iteration,
      declared_f = declared$tests$treatment$statistic,
      declared_p = declared$tests$treatment$p_value,
      classical_f = classical_f, classical_p = classical_p,
      additive_f = additive$tests$treatment$statistic,
      additive_p = additive$tests$treatment$p_value,
      stringsAsFactors = FALSE
    )
  }
  elapsed <- proc.time()[["elapsed"]] - start
  raw <- do.call(rbind, rows)
  methods <- list(
    declared_interaction = raw$declared_p,
    classical_interaction_denominator = raw$classical_p,
    additive_opt_in_misspecified = raw$additive_p
  )
  summary <- do.call(rbind, lapply(names(methods), function(method) {
    rejected <- sum(methods[[method]] <= alpha)
    uncertainty <- .binomial_summary(rejected, simulations)
    data.frame(
      method = method, simulations = simulations, rejected = rejected,
      rejection_rate = uncertainty[["rejection_rate"]],
      monte_carlo_se = uncertainty[["monte_carlo_se"]],
      confidence_low = uncertainty[["confidence_low"]],
      confidence_high = uncertainty[["confidence_high"]],
      mean_f = switch(method,
        declared_interaction = mean(raw$declared_f),
        classical_interaction_denominator = mean(raw$classical_f),
        additive_opt_in_misspecified = mean(raw$additive_f)
      ),
      stringsAsFactors = FALSE
    )
  }))
  rownames(summary) <- NULL
  out <- list(
    summary = summary,
    max_declared_classical_f_difference = max(abs(raw$declared_f - raw$classical_f)),
    max_declared_classical_p_difference = max(abs(raw$declared_p - raw$classical_p)),
    results = if (keep_results) raw else NULL,
    settings = list(
      simulations = simulations, blocks = blocks, replicates = replicates,
      interaction_sd = interaction_sd, block_sd = block_sd,
      residual_sd = residual_sd, alpha = alpha, seed = seed
    ),
    elapsed_seconds = elapsed
  )
  class(out) <- "permanova_crossed_validation"
  out
}

#' @export
print.permanova_crossed_validation <- function(x, ...) {
  cat("Fixed-by-random interaction validation\n")
  print(x$summary, row.names = FALSE)
  cat("Maximum declared/classical F difference:",
      x$max_declared_classical_f_difference, "\n")
  invisible(x)
}

.direct_blocked_permutation_reference <- function(response, data,
                                                  permutations, seed) {
  project <- function(X) {
    decomposition <- qr(X)
    Q <- qr.Q(decomposition)[, seq_len(decomposition$rank), drop = FALSE]
    list(H = tcrossprod(Q), rank = decomposition$rank)
  }
  intercept <- project(matrix(1, nrow(data), 1L))
  treatment_model <- project(model.matrix(~ treatment, data))
  additive_model <- project(model.matrix(~ treatment + block, data))
  full_model <- project(model.matrix(~ treatment * block, data))
  H_treatment <- treatment_model$H - intercept$H
  H_interaction <- full_model$H - additive_model$H
  treatment_df <- treatment_model$rank - intercept$rank
  interaction_df <- full_model$rank - additive_model$rank
  statistic_for <- function(values) {
    numerator <- sum(values * (H_treatment %*% values)) / treatment_df
    denominator <- sum(values * (H_interaction %*% values)) / interaction_df
    numerator / denominator
  }

  cell <- interaction(data$treatment, data$block, drop = TRUE)
  cell_blocks <- unname(split(seq_len(nrow(data)), cell))
  mapping <- unique(data.frame(
    cell = as.integer(cell), treatment = as.integer(data$treatment),
    block = as.integer(data$block)
  ))
  mapping <- mapping[order(mapping$block, mapping$treatment), , drop = FALSE]
  cells_by_block <- unname(split(mapping$cell, mapping$block))
  n_cells <- length(cell_blocks)
  set.seed(seed)
  orders <- matrix(rep(seq_len(n_cells), each = permutations),
                   nrow = permutations)
  for (row in seq_len(permutations)) {
    for (cells in cells_by_block) {
      orders[row, cells] <- sample(cells, length(cells))
    }
  }
  statistics <- vapply(seq_len(permutations), function(row) {
    index <- seq_len(nrow(data))
    for (target in seq_along(cell_blocks)) {
      index[cell_blocks[[target]]] <- cell_blocks[[orders[row, target]]]
    }
    statistic_for(response[index])
  }, numeric(1))
  observed <- statistic_for(response)
  list(
    statistic = observed, permutation_statistics = statistics,
    p_value = (sum(statistics >= observed - 1e-10) + 1) /
      (permutations + 1)
  )
}

#' Validate crossed randomized-block permutation inference
#'
#' Simulates a genuine treatment-by-block variance component under a null fixed
#' treatment effect. It compares the distance engine with a direct
#' response-level within-block cell permutation and the classical interaction-F
#' reference. A second scenario deliberately introduces treatment-specific
#' residual dispersion to expose exchangeability sensitivity.
#'
#' @param simulations Number of simulations per scenario.
#' @param permutations Monte Carlo permutations per fitted model.
#' @param blocks Number of random blocks.
#' @param replicates Observations per treatment-by-block cell.
#' @param interaction_sd,block_sd,residual_sd Generating standard deviations.
#' @param dispersion_ratio Impact/control residual SD ratio in the sensitivity
#'   scenario.
#' @param alpha Rejection threshold.
#' @param seed Reproducibility seed.
#' @param keep_results Retain iteration-level results.
#' @return A crossed-permutation validation report.
#' @export
simulate_crossed_permutation_validation <- function(
    simulations = 500L, permutations = 99L, blocks = 8L, replicates = 3L,
    interaction_sd = 2, block_sd = 1, residual_sd = 1,
    dispersion_ratio = 2, alpha = 0.05, seed = 1L,
    keep_results = FALSE) {
  simulations <- .positive_integer(simulations, "simulations")
  permutations <- .positive_integer(permutations, "permutations")
  blocks <- .positive_integer(blocks, "blocks")
  replicates <- .positive_integer(replicates, "replicates")
  if (blocks < 3L || replicates < 2L) {
    stop("At least three blocks and two cell replicates are required.",
         call. = FALSE)
  }
  attainable <- .safe_factorial_product(rep(2L, blocks))
  if (is.finite(attainable) && permutations >= attainable - 1L) {
    stop(
      "This validator requires a Monte Carlo sample smaller than the attainable space.",
      call. = FALSE
    )
  }
  scales <- c(interaction_sd, block_sd, residual_sd, dispersion_ratio)
  if (any(!is.finite(scales)) || any(scales < 0)) {
    stop("Generating standard deviations and dispersion ratio must be finite and non-negative.",
         call. = FALSE)
  }
  if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must lie strictly between zero and one.", call. = FALSE)
  }

  data <- expand.grid(
    replicate = factor(seq_len(replicates)),
    treatment = factor(c("control", "impact")),
    block = factor(seq_len(blocks))
  )
  cell <- interaction(data$treatment, data$block, drop = TRUE)
  scenarios <- c(
    exchangeable_null = 1,
    dispersion_heterogeneity = dispersion_ratio
  )
  raw_by_scenario <- vector("list", length(scenarios))
  names(raw_by_scenario) <- names(scenarios)
  elapsed <- numeric(length(scenarios))

  for (scenario_index in seq_along(scenarios)) {
    start <- proc.time()[["elapsed"]]
    rows <- vector("list", simulations)
    for (iteration in seq_len(simulations)) {
      iteration_seed <- as.integer(
        seed + scenario_index * 100000L + iteration
      )
      set.seed(iteration_seed)
      scale <- residual_sd * ifelse(
        data$treatment == "impact", scenarios[[scenario_index]], 1
      )
      response <- rnorm(blocks, sd = block_sd)[data$block] +
        rnorm(nlevels(cell), sd = interaction_sd)[cell] +
        rnorm(nrow(data)) * scale
      design <- compile_permanova_design(
        ~ treatment + (1 | block) + (1 | treatment:block),
        data, dist(response)
      )
      fit <- suppressWarnings(permanova_crossed(
        design, p_method = "permutation", permutations = permutations,
        seed = iteration_seed, keep_permutations = TRUE
      ))
      direct <- .direct_blocked_permutation_reference(
        response, data, permutations = permutations, seed = iteration_seed
      )
      reference <- summary(aov(response ~ treatment * block, data))[[1]]
      reference_rows <- trimws(rownames(reference))
      treatment_row <- match("treatment", reference_rows)
      interaction_row <- match("treatment:block", reference_rows)
      classical_f <- reference[["Mean Sq"]][[treatment_row]] /
        reference[["Mean Sq"]][[interaction_row]]
      classical_p <- stats::pf(
        classical_f, reference[["Df"]][[treatment_row]],
        reference[["Df"]][[interaction_row]], lower.tail = FALSE
      )
      rows[[iteration]] <- data.frame(
        iteration = iteration,
        engine_f = fit$tests$treatment$statistic,
        engine_p = fit$tests$treatment$p_value,
        direct_f = direct$statistic, direct_p = direct$p_value,
        max_permutation_statistic_difference = max(abs(
          fit$tests$treatment$permutation_statistics -
            direct$permutation_statistics
        )),
        classical_f = classical_f, classical_p = classical_p,
        stringsAsFactors = FALSE
      )
    }
    elapsed[[scenario_index]] <- proc.time()[["elapsed"]] - start
    raw_by_scenario[[scenario_index]] <- do.call(rbind, rows)
    raw_by_scenario[[scenario_index]]$scenario <- names(scenarios)[[scenario_index]]
  }
  raw <- do.call(rbind, raw_by_scenario)
  rownames(raw) <- NULL
  methods <- c("engine", "direct", "classical")
  summary <- do.call(rbind, lapply(names(scenarios), function(scenario) {
    values <- raw[raw$scenario == scenario, , drop = FALSE]
    do.call(rbind, lapply(methods, function(method) {
      p_values <- values[[paste0(method, "_p")]]
      f_values <- values[[paste0(method, "_f")]]
      rejected <- sum(p_values <= alpha)
      uncertainty <- .binomial_summary(rejected, simulations)
      data.frame(
        scenario = scenario, method = method, simulations = simulations,
        rejected = rejected,
        rejection_rate = uncertainty[["rejection_rate"]],
        monte_carlo_se = uncertainty[["monte_carlo_se"]],
        confidence_low = uncertainty[["confidence_low"]],
        confidence_high = uncertainty[["confidence_high"]],
        mean_f = mean(f_values), stringsAsFactors = FALSE
      )
    }))
  }))
  rownames(summary) <- NULL
  out <- list(
    summary = summary,
    max_engine_direct_f_difference = max(abs(raw$engine_f - raw$direct_f)),
    max_engine_direct_p_difference = max(abs(raw$engine_p - raw$direct_p)),
    max_engine_direct_permutation_statistic_difference = max(
      raw$max_permutation_statistic_difference
    ),
    max_engine_classical_f_difference = max(abs(
      raw$engine_f - raw$classical_f
    )),
    results = if (keep_results) raw else NULL,
    settings = list(
      simulations = simulations, permutations = permutations,
      blocks = blocks, replicates = replicates,
      interaction_sd = interaction_sd, block_sd = block_sd,
      residual_sd = residual_sd, dispersion_ratio = dispersion_ratio,
      alpha = alpha, seed = seed
    ),
    elapsed_seconds = setNames(elapsed, names(scenarios))
  )
  class(out) <- "permanova_crossed_permutation_validation"
  out
}

#' @export
print.permanova_crossed_permutation_validation <- function(x, ...) {
  cat("Crossed randomized-block permutation validation\n")
  print(x$summary, row.names = FALSE)
  cat("Maximum engine/direct permutation-statistic difference:",
      x$max_engine_direct_permutation_statistic_difference, "\n")
  invisible(x)
}

#' Validate balanced nested mixed-model Type I error
#'
#' Simulates treatment assigned to top-level regions with sites nested within
#' regions. The distance-based statistic is compared with ANOVA on region
#' means, and binomial uncertainty is reported for the null rejection rate.
#'
#' @param simulations Number of simulation replicates.
#' @param regions_per_treatment Top-level units per treatment.
#' @param sites_per_region Nested sites per region.
#' @param replicates Observations per site.
#' @param region_sd,site_sd,residual_sd Generating standard deviations.
#' @param alpha Rejection threshold.
#' @param seed Reproducibility seed.
#' @param keep_results Retain iteration-level results.
#' @return A nested-validation report.
#' @export
simulate_nested_validation <- function(
    simulations = 500L, regions_per_treatment = 4L,
    sites_per_region = 2L, replicates = 3L,
    region_sd = 1.5, site_sd = 1, residual_sd = 0.5,
    alpha = 0.05, seed = 1L, keep_results = FALSE) {
  simulations <- .positive_integer(simulations, "simulations")
  regions_per_treatment <- .positive_integer(
    regions_per_treatment, "regions_per_treatment"
  )
  sites_per_region <- .positive_integer(sites_per_region, "sites_per_region")
  replicates <- .positive_integer(replicates, "replicates")
  if (regions_per_treatment < 2L || replicates < 2L) {
    stop("At least two regions per treatment and two observations per site are required.", call. = FALSE)
  }
  scales <- c(region_sd, site_sd, residual_sd)
  if (any(!is.finite(scales)) || any(scales < 0)) {
    stop("Generating standard deviations must be finite and non-negative.", call. = FALSE)
  }
  if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must lie strictly between zero and one.", call. = FALSE)
  }

  regions <- 2L * regions_per_treatment
  data <- expand.grid(
    replicate = factor(seq_len(replicates)),
    site = factor(seq_len(sites_per_region)),
    region = factor(seq_len(regions))
  )
  data$treatment <- factor(
    ifelse(as.integer(data$region) <= regions_per_treatment, "control", "impact")
  )
  rows <- vector("list", simulations)
  start <- proc.time()[["elapsed"]]
  for (iteration in seq_len(simulations)) {
    set.seed(as.integer(seed + iteration))
    region_effect <- rnorm(regions, sd = region_sd)
    nested_site <- interaction(data$region, data$site, drop = TRUE)
    site_effect <- rnorm(nlevels(nested_site), sd = site_sd)
    response <- region_effect[data$region] + site_effect[nested_site] +
      rnorm(nrow(data), sd = residual_sd)
    design <- compile_permanova_design(
      ~ treatment + (1 | region/site), data, dist(response)
    )
    fit <- permanova_nested(design, permutations = 0)
    region_means <- aggregate(
      response, list(region = data$region, treatment = data$treatment), mean
    )
    names(region_means)[[3L]] <- "response"
    reference <- summary(aov(response ~ treatment, region_means))[[1]]
    classical_f <- unname(reference[["F value"]][[1L]])
    classical_p <- unname(reference[["Pr(>F)"]][[1L]])
    rows[[iteration]] <- data.frame(
      iteration = iteration,
      nested_f = fit$tests$treatment$statistic,
      nested_reference_p = stats::pf(
        fit$tests$treatment$statistic,
        fit$tests$treatment$numerator$df,
        fit$tests$treatment$denominator$df,
        lower.tail = FALSE
      ),
      classical_f = classical_f, classical_p = classical_p,
      stringsAsFactors = FALSE
    )
  }
  elapsed <- proc.time()[["elapsed"]] - start
  raw <- do.call(rbind, rows)
  methods <- list(
    nested_exact_denominator = raw$nested_reference_p,
    classical_region_mean_anova = raw$classical_p
  )
  summary <- do.call(rbind, lapply(names(methods), function(method) {
    rejected <- sum(methods[[method]] <= alpha)
    uncertainty <- .binomial_summary(rejected, simulations)
    data.frame(
      method = method, simulations = simulations, rejected = rejected,
      rejection_rate = uncertainty[["rejection_rate"]],
      monte_carlo_se = uncertainty[["monte_carlo_se"]],
      confidence_low = uncertainty[["confidence_low"]],
      confidence_high = uncertainty[["confidence_high"]],
      mean_f = if (method == "nested_exact_denominator") {
        mean(raw$nested_f)
      } else mean(raw$classical_f),
      stringsAsFactors = FALSE
    )
  }))
  rownames(summary) <- NULL
  out <- list(
    summary = summary,
    max_nested_classical_f_difference = max(abs(raw$nested_f - raw$classical_f)),
    max_nested_classical_p_difference = max(abs(raw$nested_reference_p - raw$classical_p)),
    results = if (keep_results) raw else NULL,
    settings = list(
      simulations = simulations, regions_per_treatment = regions_per_treatment,
      sites_per_region = sites_per_region, replicates = replicates,
      region_sd = region_sd, site_sd = site_sd, residual_sd = residual_sd,
      alpha = alpha, seed = seed
    ),
    elapsed_seconds = elapsed
  )
  class(out) <- "permanova_nested_validation"
  out
}

#' @export
print.permanova_nested_validation <- function(x, ...) {
  cat("Balanced nested mixed-model validation\n")
  print(x$summary, row.names = FALSE)
  cat("Maximum nested/classical F difference:",
      x$max_nested_classical_f_difference, "\n")
  invisible(x)
}

#' Validate level-specific balanced nested inference
#'
#' Simulates null fixed effects assigned to regions, sites within regions, and
#' observations within sites. Each distance-based pseudo-F and its reference F
#' tail probability are compared with the corresponding classical nested
#' ANOVA calculation.
#'
#' @param simulations Number of simulation replicates.
#' @param regions_per_treatment Regions per top-level treatment.
#' @param replicates Observations per site.
#' @param region_sd,site_sd,residual_sd Generating standard deviations.
#' @param alpha Rejection threshold.
#' @param seed Reproducibility seed.
#' @param keep_results Retain iteration-level results.
#' @return A level-specific nested-validation report.
#' @export
simulate_nested_level_validation <- function(
    simulations = 500L, regions_per_treatment = 4L, replicates = 3L,
    region_sd = 1.5, site_sd = 1, residual_sd = 0.5,
    alpha = 0.05, seed = 1L, keep_results = FALSE) {
  simulations <- .positive_integer(simulations, "simulations")
  regions_per_treatment <- .positive_integer(
    regions_per_treatment, "regions_per_treatment"
  )
  replicates <- .positive_integer(replicates, "replicates")
  if (regions_per_treatment < 2L || replicates < 2L) {
    stop("At least two regions per treatment and two observations per site are required.",
         call. = FALSE)
  }
  scales <- c(region_sd, site_sd, residual_sd)
  if (any(!is.finite(scales)) || any(scales < 0)) {
    stop("Generating standard deviations must be finite and non-negative.",
         call. = FALSE)
  }
  if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must lie strictly between zero and one.", call. = FALSE)
  }

  regions <- 2L * regions_per_treatment
  data <- expand.grid(
    replicate = factor(seq_len(replicates)), site = factor(1:2),
    region = factor(seq_len(regions))
  )
  data$region_treatment <- factor(
    ifelse(as.integer(data$region) <= regions_per_treatment,
           "control", "impact")
  )
  data$site_treatment <- data$site
  data$observation_treatment <- data$replicate
  nested_site <- interaction(data$region, data$site, drop = TRUE)
  rows <- vector("list", simulations)
  start <- proc.time()[["elapsed"]]

  for (iteration in seq_len(simulations)) {
    set.seed(as.integer(seed + iteration))
    response <- rnorm(regions, sd = region_sd)[data$region] +
      rnorm(nlevels(nested_site), sd = site_sd)[nested_site] +
      rnorm(nrow(data), sd = residual_sd)
    design <- compile_permanova_design(
      ~ region_treatment + site_treatment + observation_treatment +
        (1 | region/site),
      data, dist(response)
    )
    fit <- permanova_nested(design, permutations = 0)

    region_means <- aggregate(
      response,
      list(region = data$region, region_treatment = data$region_treatment),
      mean
    )
    names(region_means)[[3L]] <- "response"
    site_means <- aggregate(
      response,
      list(
        region = data$region, site = data$site,
        site_treatment = data$site_treatment
      ),
      mean
    )
    names(site_means)[[4L]] <- "response"
    region_reference <- summary(aov(
      response ~ region_treatment, region_means
    ))[[1]]
    site_reference <- summary(aov(
      response ~ region + site_treatment, site_means
    ))[[1]]
    observation_reference <- summary(aov(
      response ~ nested_site + observation_treatment, data
    ))[[1]]
    references <- list(
      region = region_reference["region_treatment", , drop = FALSE],
      site = site_reference["site_treatment", , drop = FALSE],
      observation = observation_reference[
        "observation_treatment", , drop = FALSE
      ]
    )
    terms <- c(
      region = "region_treatment", site = "site_treatment",
      observation = "observation_treatment"
    )
    values <- lapply(names(terms), function(level) {
      specification <- fit$tests[[terms[[level]]]]
      reference <- references[[level]]
      c(
        nested_f = specification$statistic,
        nested_p = stats::pf(
          specification$statistic, specification$numerator$df,
          specification$denominator$df, lower.tail = FALSE
        ),
        classical_f = unname(reference[["F value"]][[1L]]),
        classical_p = unname(reference[["Pr(>F)"]][[1L]])
      )
    })
    names(values) <- names(terms)
    rows[[iteration]] <- data.frame(
      iteration = iteration,
      region_nested_f = values$region[["nested_f"]],
      region_nested_p = values$region[["nested_p"]],
      region_classical_f = values$region[["classical_f"]],
      region_classical_p = values$region[["classical_p"]],
      site_nested_f = values$site[["nested_f"]],
      site_nested_p = values$site[["nested_p"]],
      site_classical_f = values$site[["classical_f"]],
      site_classical_p = values$site[["classical_p"]],
      observation_nested_f = values$observation[["nested_f"]],
      observation_nested_p = values$observation[["nested_p"]],
      observation_classical_f = values$observation[["classical_f"]],
      observation_classical_p = values$observation[["classical_p"]],
      stringsAsFactors = FALSE
    )
  }

  elapsed <- proc.time()[["elapsed"]] - start
  raw <- do.call(rbind, rows)
  levels <- c("region", "site", "observation")
  summary <- do.call(rbind, lapply(levels, function(level) {
    do.call(rbind, lapply(c("nested", "classical"), function(method) {
      p_values <- raw[[paste0(level, "_", method, "_p")]]
      f_values <- raw[[paste0(level, "_", method, "_f")]]
      rejected <- sum(p_values <= alpha)
      uncertainty <- .binomial_summary(rejected, simulations)
      data.frame(
        level = level, method = method, simulations = simulations,
        rejected = rejected,
        rejection_rate = uncertainty[["rejection_rate"]],
        monte_carlo_se = uncertainty[["monte_carlo_se"]],
        confidence_low = uncertainty[["confidence_low"]],
        confidence_high = uncertainty[["confidence_high"]],
        mean_f = mean(f_values), stringsAsFactors = FALSE
      )
    }))
  }))
  rownames(summary) <- NULL
  max_f_difference <- setNames(vapply(levels, function(level) {
    max(abs(
      raw[[paste0(level, "_nested_f")]] -
        raw[[paste0(level, "_classical_f")]]
    ))
  }, numeric(1)), levels)
  max_p_difference <- setNames(vapply(levels, function(level) {
    max(abs(
      raw[[paste0(level, "_nested_p")]] -
        raw[[paste0(level, "_classical_p")]]
    ))
  }, numeric(1)), levels)
  out <- list(
    summary = summary, max_f_difference = max_f_difference,
    max_p_difference = max_p_difference,
    results = if (keep_results) raw else NULL,
    settings = list(
      simulations = simulations,
      regions_per_treatment = regions_per_treatment,
      sites_per_region = 2L, replicates = replicates,
      region_sd = region_sd, site_sd = site_sd, residual_sd = residual_sd,
      alpha = alpha, seed = seed
    ),
    elapsed_seconds = elapsed
  )
  class(out) <- "permanova_nested_level_validation"
  out
}

#' @export
print.permanova_nested_level_validation <- function(x, ...) {
  cat("Level-specific balanced nested validation\n")
  print(x$summary, row.names = FALSE)
  cat("Maximum nested/classical F differences:\n")
  print(x$max_f_difference)
  invisible(x)
}

.dispersion_simulation_case <- function(structure, scenario, seed) {
  set.seed(seed)
  dimensions <- 4L
  dispersion_ratio <- if (scenario == "pure_dispersion") 2 else 1
  location_shift <- if (scenario == "location_only") 1.5 else 0

  if (structure == "one_way") {
    data <- data.frame(treatment = factor(rep(c("control", "impact"), each = 12)))
    scale <- ifelse(data$treatment == "impact", dispersion_ratio, 1)
    response <- matrix(rnorm(nrow(data) * dimensions), ncol = dimensions) * scale
    response[, 1L] <- response[, 1L] +
      location_shift * (data$treatment == "impact")
    formula <- ~ treatment
  } else if (structure == "random_intercept") {
    data <- expand.grid(replicate = factor(1:3), cluster = factor(1:8))
    data$treatment <- factor(
      ifelse(as.integer(data$cluster) <= 4, "control", "impact")
    )
    cluster_effect <- matrix(rnorm(8 * dimensions, sd = 0.7), 8, dimensions)
    scale <- ifelse(data$treatment == "impact", dispersion_ratio, 1)
    response <- cluster_effect[data$cluster, ] +
      matrix(rnorm(nrow(data) * dimensions), ncol = dimensions) * scale
    response[, 1L] <- response[, 1L] +
      location_shift * (data$treatment == "impact")
    formula <- ~ treatment + (1 | cluster)
  } else if (structure == "nested") {
    data <- expand.grid(
      replicate = factor(1:3), site = factor(1:2), region = factor(1:6)
    )
    data$treatment <- data$site
    nested_site <- interaction(data$region, data$site, drop = TRUE)
    region_effect <- matrix(rnorm(6 * dimensions, sd = 0.7), 6, dimensions)
    site_effect <- matrix(rnorm(12 * dimensions, sd = 0.45), 12, dimensions)
    scale <- ifelse(data$treatment == "2", dispersion_ratio, 1)
    response <- region_effect[data$region, ] + site_effect[nested_site, ] +
      matrix(rnorm(nrow(data) * dimensions), ncol = dimensions) * scale
    response[, 1L] <- response[, 1L] +
      location_shift * (data$treatment == "2")
    formula <- ~ treatment + (1 | region/site)
  } else if (structure == "crossed") {
    data <- expand.grid(
      replicate = factor(1:3), treatment = factor(c("control", "impact")),
      block = factor(1:8)
    )
    cell <- interaction(data$treatment, data$block, drop = TRUE)
    block_effect <- matrix(rnorm(8 * dimensions, sd = 0.7), 8, dimensions)
    cell_effect <- matrix(rnorm(16 * dimensions, sd = 0.55), 16, dimensions)
    scale <- ifelse(data$treatment == "impact", dispersion_ratio, 1)
    response <- block_effect[data$block, ] + cell_effect[cell, ] +
      matrix(rnorm(nrow(data) * dimensions), ncol = dimensions) * scale
    response[, 1L] <- response[, 1L] +
      location_shift * (data$treatment == "impact")
    formula <- ~ treatment + (1 | block) + (1 | treatment:block)
  } else {
    stop("Unknown dispersion-validation structure.", call. = FALSE)
  }
  list(data = data, response = response, formula = formula)
}

.dispersion_classical_f <- function(structure, distances, data) {
  if (structure == "one_way") {
    table <- summary(aov(distances ~ treatment, data))[[1]]
    return(unname(table[1L, "F value"]))
  }
  if (structure == "random_intercept") {
    means <- aggregate(
      distances,
      list(cluster = data$cluster, treatment = data$treatment), mean
    )
    table <- summary(aov(x ~ treatment, means))[[1]]
    return(unname(table[1L, "F value"]))
  }
  if (structure == "nested") {
    means <- aggregate(
      distances,
      list(region = data$region, site = data$site,
           treatment = data$treatment), mean
    )
    table <- summary(aov(x ~ region + treatment, means))[[1]]
    return(unname(table["treatment", "F value"]))
  }
  table <- summary(aov(
    distances ~ treatment * block, data
  ))[[1]]
  rows <- trimws(rownames(table))
  unname(table[match("treatment", rows), "Mean Sq"] /
           table[match("treatment:block", rows), "Mean Sq"])
}

#' Validate design-aware multivariate dispersion inference
#'
#' Evaluates equal-dispersion null calibration, invariance to a pure location
#' shift, and power against a dispersion ratio across one-way, random-intercept,
#' nested, and canonical randomized-block designs. Every observed statistic is
#' compared with an independently fitted scalar classical reference.
#'
#' @param simulations Number of simulations per structure and scenario.
#' @param permutations Monte Carlo residual permutations per fit.
#' @param alpha Rejection threshold.
#' @param seed Reproducibility seed.
#' @param type Group spatial center: `"median"` or `"centroid"`.
#' @param keep_results Retain iteration-level results.
#' @return A dispersion-validation report.
#' @export
simulate_dispersion_validation <- function(
    simulations = 200L, permutations = 99L, alpha = 0.05,
    seed = 20260718L, type = c("median", "centroid"),
    keep_results = FALSE) {
  type <- match.arg(type)
  simulations <- .positive_integer(simulations, "simulations")
  permutations <- .positive_integer(permutations, "permutations")
  if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must lie strictly between zero and one.", call. = FALSE)
  }
  structures <- c("one_way", "random_intercept", "nested", "crossed")
  scenarios <- c("exchangeable_null", "location_only", "pure_dispersion")
  rows <- vector("list", length(structures) * length(scenarios) * simulations)
  cursor <- 0L
  start <- proc.time()[["elapsed"]]

  for (structure_index in seq_along(structures)) {
    structure <- structures[[structure_index]]
    for (scenario_index in seq_along(scenarios)) {
      scenario <- scenarios[[scenario_index]]
      generator_scenario <- switch(
        scenario, exchangeable_null = "null",
        location_only = "location_only", pure_dispersion = "pure_dispersion"
      )
      for (iteration in seq_len(simulations)) {
        iteration_seed <- as.integer(
          seed + structure_index * 1000000L + scenario_index * 10000L + iteration
        )
        generated <- .dispersion_simulation_case(
          structure, generator_scenario, iteration_seed
        )
        design <- compile_permanova_design(
          generated$formula, generated$data, dist(generated$response)
        )
        fit <- suppressWarnings(permanova_dispersion(
          design, "treatment", type = type,
          permutations = permutations, seed = iteration_seed + 5000L
        ))
        reference_f <- .dispersion_classical_f(
          structure, fit$distances, generated$data
        )
        cursor <- cursor + 1L
        rows[[cursor]] <- data.frame(
          structure = structure, scenario = scenario, iteration = iteration,
          engine_f = fit$test$statistic,
          reference_f = reference_f,
          p_value = fit$test$p_value,
          rejected = is.finite(fit$test$p_value) && fit$test$p_value <= alpha,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  raw <- do.call(rbind, rows)
  summary <- do.call(rbind, lapply(structures, function(structure) {
    do.call(rbind, lapply(scenarios, function(scenario) {
      selected <- raw$structure == structure & raw$scenario == scenario
      rejected <- sum(raw$rejected[selected])
      uncertainty <- .binomial_summary(rejected, simulations)
      data.frame(
        structure = structure, scenario = scenario,
        estimand = if (scenario == "pure_dispersion") "power" else "type_I_error",
        simulations = simulations, rejected = rejected,
        rejection_rate = uncertainty[["rejection_rate"]],
        monte_carlo_se = uncertainty[["monte_carlo_se"]],
        confidence_low = uncertainty[["confidence_low"]],
        confidence_high = uncertainty[["confidence_high"]],
        mean_f = mean(raw$engine_f[selected]), stringsAsFactors = FALSE
      )
    }))
  }))
  rownames(summary) <- NULL
  max_difference <- max(abs(raw$engine_f - raw$reference_f))
  out <- list(
    summary = summary,
    max_engine_classical_f_difference = max_difference,
    results = if (keep_results) raw else NULL,
    settings = list(
      simulations = simulations, permutations = permutations,
      structures = structures, scenarios = scenarios,
      dispersion_ratio = 2, location_shift = 1.5,
      center = type, alpha = alpha, seed = seed
    ),
    elapsed_seconds = proc.time()[["elapsed"]] - start
  )
  class(out) <- "permanova_dispersion_validation"
  out
}

#' @export
print.permanova_dispersion_validation <- function(x, ...) {
  cat("Design-aware dispersion validation\n")
  print(x$summary, row.names = FALSE)
  cat("Maximum engine/classical F difference:",
      format(x$max_engine_classical_f_difference), "\n")
  invisible(x)
}

#' Validate balanced crossed-nested hybrid inference
#'
#' Simulates a fixed treatment crossed with regions and sites nested within
#' regions. The engine is compared with ANOVA on treatment-by-region profile
#' means and an independent response-level implementation of every restricted
#' treatment-profile permutation.
#'
#' @param simulations Number of simulations per null and location scenario.
#' @param permutations Monte Carlo treatment-profile permutations per fit.
#' @param regions,sites,replicates Balanced hierarchy dimensions.
#' @param region_sd,site_sd,interaction_sd,nested_interaction_sd,residual_sd
#'   Generating standard deviations.
#' @param effect Fixed treatment shift in the power scenario.
#' @param alpha Rejection threshold.
#' @param seed Reproducibility seed.
#' @param keep_results Retain iteration-level results.
#' @return A hybrid-validation report.
#' @export
simulate_hybrid_validation <- function(
    simulations = 500L, permutations = 99L, regions = 8L,
    sites = 3L, replicates = 2L, region_sd = 1,
    site_sd = 0.7, interaction_sd = 1.2,
    nested_interaction_sd = 0.6, residual_sd = 0.4,
    effect = 1, alpha = 0.05, seed = 20260718L,
    keep_results = FALSE) {
  simulations <- .positive_integer(simulations, "simulations")
  permutations <- .positive_integer(permutations, "permutations")
  regions <- .positive_integer(regions, "regions")
  sites <- .positive_integer(sites, "sites")
  replicates <- .positive_integer(replicates, "replicates")
  if (regions < 2L || sites < 2L || replicates < 2L) {
    stop("Hybrid validation requires at least two regions, sites, and replicates.",
         call. = FALSE)
  }
  scales <- c(
    region_sd, site_sd, interaction_sd, nested_interaction_sd, residual_sd
  )
  if (any(!is.finite(scales)) || any(scales < 0) ||
      !is.finite(effect)) {
    stop("Generating effects must be finite and standard deviations non-negative.",
         call. = FALSE)
  }
  if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must lie strictly between zero and one.", call. = FALSE)
  }
  data <- expand.grid(
    replicate = factor(seq_len(replicates)),
    treatment = factor(c("control", "impact")),
    site = factor(seq_len(sites)), region = factor(seq_len(regions))
  )
  region_site <- interaction(data$region, data$site, drop = TRUE)
  treatment_region <- interaction(
    data$treatment, data$region, drop = TRUE
  )
  treatment_site <- interaction(
    data$treatment, data$region, data$site, drop = TRUE
  )
  formula <- ~ treatment + (1 | region/site) +
    (1 | treatment:region) + (1 | treatment:region:site)
  scenarios <- c(exchangeable_null = 0, location_shift = effect)
  rows <- vector("list", length(scenarios) * simulations)
  cursor <- 0L
  start <- proc.time()[["elapsed"]]

  for (scenario_index in seq_along(scenarios)) {
    scenario <- names(scenarios)[[scenario_index]]
    fixed_effect <- scenarios[[scenario_index]]
    for (iteration in seq_len(simulations)) {
      iteration_seed <- as.integer(
        seed + scenario_index * 100000L + iteration
      )
      set.seed(iteration_seed)
      response <-
        rnorm(regions, sd = region_sd)[data$region] +
        rnorm(regions * sites, sd = site_sd)[region_site] +
        rnorm(2L * regions, sd = interaction_sd)[treatment_region] +
        rnorm(2L * regions * sites, sd = nested_interaction_sd)[treatment_site] +
        rnorm(nrow(data), sd = residual_sd) +
        fixed_effect * (data$treatment == "impact")
      design <- compile_permanova_design(formula, data, dist(response))
      fit <- suppressWarnings(permanova_hybrid(
        design, permutations = permutations,
        seed = iteration_seed + 5000L, keep_permutations = TRUE
      ))
      structure <- .validate_hybrid_design(design)
      set.seed(iteration_seed + 5000L)
      plan <- .hybrid_permutation_plan(
        structure, permutations = permutations, max_exact = 100000L
      )
      profile_f <- function(values) {
        means <- aggregate(
          values,
          list(treatment = data$treatment, region = data$region), mean
        )
        table <- summary(aov(x ~ treatment * region, means))[[1]]
        table_rows <- trimws(rownames(table))
        unname(table[match("treatment", table_rows), "Mean Sq"] /
                 table[match("treatment:region", table_rows), "Mean Sq"])
      }
      classical_f <- profile_f(response)
      direct_statistics <- vapply(
        seq_len(nrow(plan$orders)), function(row) {
          index <- .block_permutation_index(
            structure$profile_blocks, plan$orders[row, ]
          )
          profile_f(response[index])
        }, numeric(1)
      )
      extreme <- sum(direct_statistics >= classical_f - 1e-10)
      direct_p <- if (plan$exact) {
        extreme / length(direct_statistics)
      } else {
        (extreme + 1) / (length(direct_statistics) + 1)
      }
      classical_p <- stats::pf(
        classical_f, nlevels(data$treatment) - 1L,
        (nlevels(data$treatment) - 1L) * (regions - 1L),
        lower.tail = FALSE
      )
      cursor <- cursor + 1L
      rows[[cursor]] <- data.frame(
        scenario = scenario, iteration = iteration,
        engine_f = fit$test$statistic, engine_p = fit$test$p_value,
        direct_f = classical_f, direct_p = direct_p,
        classical_f = classical_f, classical_p = classical_p,
        max_permutation_difference = max(abs(
          fit$test$permutation_statistics - direct_statistics
        )), stringsAsFactors = FALSE
      )
    }
  }
  raw <- do.call(rbind, rows)
  methods <- c("engine", "direct", "classical")
  summary <- do.call(rbind, lapply(names(scenarios), function(scenario) {
    selected <- raw$scenario == scenario
    do.call(rbind, lapply(methods, function(method) {
      p_values <- raw[[paste0(method, "_p")]][selected]
      rejected <- sum(p_values <= alpha)
      uncertainty <- .binomial_summary(rejected, simulations)
      data.frame(
        scenario = scenario, method = method,
        estimand = if (scenario == "exchangeable_null") "type_I_error" else "power",
        simulations = simulations, rejected = rejected,
        rejection_rate = uncertainty[["rejection_rate"]],
        monte_carlo_se = uncertainty[["monte_carlo_se"]],
        confidence_low = uncertainty[["confidence_low"]],
        confidence_high = uncertainty[["confidence_high"]],
        stringsAsFactors = FALSE
      )
    }))
  }))
  rownames(summary) <- NULL
  out <- list(
    summary = summary,
    max_engine_direct_f_difference = max(abs(raw$engine_f - raw$direct_f)),
    max_engine_direct_p_difference = max(abs(raw$engine_p - raw$direct_p)),
    max_engine_direct_permutation_statistic_difference =
      max(raw$max_permutation_difference),
    max_engine_classical_f_difference =
      max(abs(raw$engine_f - raw$classical_f)),
    results = if (keep_results) raw else NULL,
    settings = list(
      simulations = simulations, permutations = permutations,
      regions = regions, sites = sites, replicates = replicates,
      region_sd = region_sd, site_sd = site_sd,
      interaction_sd = interaction_sd,
      nested_interaction_sd = nested_interaction_sd,
      residual_sd = residual_sd, effect = effect,
      alpha = alpha, seed = seed
    ),
    elapsed_seconds = proc.time()[["elapsed"]] - start
  )
  class(out) <- "permanova_hybrid_validation"
  out
}

#' @export
print.permanova_hybrid_validation <- function(x, ...) {
  cat("Balanced crossed-nested hybrid validation\n")
  print(x$summary, row.names = FALSE)
  cat("Maximum engine/direct F difference:",
      format(x$max_engine_direct_f_difference), "\n")
  cat("Maximum engine/direct permutation-statistic difference:",
      format(x$max_engine_direct_permutation_statistic_difference), "\n")
  invisible(x)
}

#' Validate unbalanced fixed-effect inference
#'
#' Simulates a fixed unequal-cell additive design and compares the Type II
#' distance engine with nested scalar linear models and an independent
#' response-level Freedman--Lane implementation under the null and a location
#' alternative.
#'
#' @param simulations Number of simulations per scenario.
#' @param permutations Monte Carlo residual permutations per fit.
#' @param effect Fixed A-factor effect in the power scenario.
#' @param residual_sd Generating residual standard deviation.
#' @param alpha Rejection threshold.
#' @param seed Reproducibility seed.
#' @param keep_results Retain iteration-level results.
#' @return An unbalanced-validation report.
#' @export
simulate_unbalanced_validation <- function(
    simulations = 500L, permutations = 99L, effect = 0.8,
    residual_sd = 1, alpha = 0.05, seed = 20260718L,
    keep_results = FALSE) {
  simulations <- .positive_integer(simulations, "simulations")
  permutations <- .positive_integer(permutations, "permutations")
  if (!is.finite(effect) || !is.finite(residual_sd) || residual_sd <= 0) {
    stop("`effect` must be finite and `residual_sd` positive.", call. = FALSE)
  }
  if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must lie strictly between zero and one.", call. = FALSE)
  }
  cells <- expand.grid(
    A = factor(c("a1", "a2")), B = factor(c("b1", "b2", "b3"))
  )
  cell_counts <- c(3L, 5L, 4L, 2L, 6L, 3L)
  data <- cells[rep(seq_len(nrow(cells)), cell_counts), , drop = FALSE]
  rownames(data) <- NULL
  B_effect <- c(b1 = -0.7, b2 = 0.1, b3 = 0.8)
  scenarios <- c(exchangeable_null = 0, location_shift = effect)
  rows <- vector("list", length(scenarios) * simulations)
  cursor <- 0L
  start <- proc.time()[["elapsed"]]

  for (scenario_index in seq_along(scenarios)) {
    scenario <- names(scenarios)[[scenario_index]]
    A_effect <- scenarios[[scenario_index]]
    for (iteration in seq_len(simulations)) {
      iteration_seed <- as.integer(
        seed + scenario_index * 100000L + iteration
      )
      set.seed(iteration_seed)
      response <- A_effect * (data$A == "a2") + B_effect[data$B] +
        rnorm(nrow(data), sd = residual_sd)
      additive_formula <- stats::reformulate(c("B", "A"))
      reduced_formula <- stats::reformulate("B")
      design <- compile_permanova_design(
        additive_formula, data, dist(response)
      )
      fit <- permanova_unbalanced(
        design, hypothesis = "II", terms = "A",
        permutations = permutations, seed = iteration_seed + 5000L,
        keep_permutations = TRUE
      )
      full <- stats::lm(
        stats::reformulate(c("B", "A"), response = "response"), data
      )
      reduced <- stats::lm(
        stats::reformulate("B", response = "response"), data
      )
      full_rss <- sum(stats::residuals(full) ^ 2)
      reduced_rss <- sum(stats::residuals(reduced) ^ 2)
      numerator_df <- stats::df.residual(reduced) - stats::df.residual(full)
      classical_f <- ((reduced_rss - full_rss) / numerator_df) /
        (full_rss / stats::df.residual(full))
      classical_p <- stats::pf(
        classical_f, numerator_df, stats::df.residual(full),
        lower.tail = FALSE
      )

      set.seed(iteration_seed + 5000L)
      plan <- .permutation_orders(
        nrow(data), permutations = permutations, max_exact = 100000L
      )
      fitted_reduced <- stats::fitted(reduced)
      residual_reduced <- stats::residuals(reduced)
      direct_statistics <- apply(plan$orders, 1L, function(index) {
        permuted_response <- fitted_reduced + residual_reduced[index]
        permutation_data <- data
        permutation_data$permuted_response <- permuted_response
        permuted_full <- stats::lm(
          stats::reformulate(
            c("B", "A"), response = "permuted_response"
          ),
          permutation_data
        )
        permuted_reduced <- stats::lm(
          stats::reformulate("B", response = "permuted_response"),
          permutation_data
        )
        rss_full <- sum(stats::residuals(permuted_full) ^ 2)
        rss_reduced <- sum(stats::residuals(permuted_reduced) ^ 2)
        ((rss_reduced - rss_full) / numerator_df) /
          (rss_full / stats::df.residual(permuted_full))
      })
      extreme <- sum(direct_statistics >= classical_f - 1e-10)
      direct_p <- if (plan$exact) {
        extreme / length(direct_statistics)
      } else {
        (extreme + 1) / (length(direct_statistics) + 1)
      }
      cursor <- cursor + 1L
      rows[[cursor]] <- data.frame(
        scenario = scenario, iteration = iteration,
        engine_f = fit$tests$A$statistic,
        engine_p = fit$tests$A$p_value,
        direct_f = classical_f, direct_p = direct_p,
        classical_f = classical_f, classical_p = classical_p,
        max_permutation_difference = max(abs(
          fit$tests$A$permutation_statistics - direct_statistics
        )), stringsAsFactors = FALSE
      )
    }
  }
  raw <- do.call(rbind, rows)
  methods <- c("engine", "direct", "classical")
  summary <- do.call(rbind, lapply(names(scenarios), function(scenario) {
    selected <- raw$scenario == scenario
    do.call(rbind, lapply(methods, function(method) {
      p_values <- raw[[paste0(method, "_p")]][selected]
      rejected <- sum(p_values <= alpha)
      uncertainty <- .binomial_summary(rejected, simulations)
      data.frame(
        scenario = scenario, method = method,
        estimand = if (scenario == "exchangeable_null") "type_I_error" else "power",
        simulations = simulations, rejected = rejected,
        rejection_rate = uncertainty[["rejection_rate"]],
        monte_carlo_se = uncertainty[["monte_carlo_se"]],
        confidence_low = uncertainty[["confidence_low"]],
        confidence_high = uncertainty[["confidence_high"]],
        stringsAsFactors = FALSE
      )
    }))
  }))
  rownames(summary) <- NULL
  out <- list(
    summary = summary,
    max_engine_direct_f_difference = max(abs(raw$engine_f - raw$direct_f)),
    max_engine_direct_p_difference = max(abs(raw$engine_p - raw$direct_p)),
    max_engine_direct_permutation_statistic_difference =
      max(raw$max_permutation_difference),
    max_engine_classical_f_difference =
      max(abs(raw$engine_f - raw$classical_f)),
    results = if (keep_results) raw else NULL,
    settings = list(
      simulations = simulations, permutations = permutations,
      cell_counts = cell_counts, effect = effect,
      residual_sd = residual_sd, alpha = alpha, seed = seed,
      hypothesis = "II", tested_term = "A", formula = "~ B + A"
    ),
    elapsed_seconds = proc.time()[["elapsed"]] - start
  )
  class(out) <- "permanova_unbalanced_validation"
  out
}

#' @export
print.permanova_unbalanced_validation <- function(x, ...) {
  cat("Unbalanced fixed-effect validation\n")
  print(x$summary, row.names = FALSE)
  cat("Maximum engine/direct F difference:",
      format(x$max_engine_direct_f_difference), "\n")
  cat("Maximum engine/direct permutation-statistic difference:",
      format(x$max_engine_direct_permutation_statistic_difference), "\n")
  invisible(x)
}

#' Validate unequal-replication random-intercept inference
#'
#' Compares the Phase 11 distance engine with an independent scalar
#' cluster-centroid Freedman--Lane implementation under a null and location
#' alternative. Cluster sizes are 2, 3, and 5, with both treatments represented
#' within every size stratum.
#'
#' @param simulations Number of simulations per scenario.
#' @param permutations Monte Carlo permutations per fit.
#' @param effect Cluster-level treatment shift in the power scenario.
#' @param cluster_sd,residual_sd Generating standard deviations.
#' @param alpha Rejection threshold.
#' @param seed Reproducibility seed.
#' @param keep_results Retain iteration-level results.
#' @return A Phase 11 validation report.
#' @export
simulate_unequal_random_intercept_validation <- function(
    simulations = 500L, permutations = 99L, effect = 1,
    cluster_sd = 1, residual_sd = 1, alpha = 0.05,
    seed = 20260719L, keep_results = FALSE) {
  simulations <- .positive_integer(simulations, "simulations")
  permutations <- .positive_integer(permutations, "permutations")
  if (any(!is.finite(c(effect, cluster_sd, residual_sd))) ||
      cluster_sd < 0 || residual_sd <= 0) {
    stop("Simulation effects must be finite and standard deviations valid.",
         call. = FALSE)
  }
  if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must lie strictly between zero and one.", call. = FALSE)
  }
  sizes <- rep(c(2L, 3L, 5L), each = 4L)
  unit_treatment <- rep(rep(c("control", "control", "impact", "impact")), 3L)
  data <- data.frame(
    site = factor(rep(seq_along(sizes), sizes)),
    treatment = factor(rep(unit_treatment, sizes),
                       levels = c("control", "impact"))
  )
  scenarios <- c(exchangeable_null = 0, location_shift = effect)
  rows <- vector("list", length(scenarios) * simulations)
  cursor <- 0L
  start <- proc.time()[["elapsed"]]
  for (scenario_index in seq_along(scenarios)) {
    scenario <- names(scenarios)[[scenario_index]]
    shift <- scenarios[[scenario_index]]
    for (iteration in seq_len(simulations)) {
      iteration_seed <- as.integer(seed + scenario_index * 100000L + iteration)
      set.seed(iteration_seed)
      unit_effect <- rnorm(length(sizes), sd = cluster_sd)
      response <- shift * (data$treatment == "impact") +
        unit_effect[data$site] + rnorm(nrow(data), sd = residual_sd)
      design <- compile_permanova_design(
        ~ treatment + (1 | site), data, dist(response)
      )
      fit <- permanova_unequal_random_intercept(
        design, permutations = permutations, seed = iteration_seed + 5000L,
        keep_permutations = TRUE
      )
      means <- as.numeric(aggregate(response, list(data$site), mean)$x)
      classical_table <- summary(aov(means ~ factor(unit_treatment)))[[1]]
      direct_f <- unname(classical_table[1L, "F value"])

      structure <- .validate_unequal_random_intercept(design)
      set.seed(iteration_seed + 5000L)
      plan <- .stratified_permutation_plan(
        structure$size_strata, permutations, 100000L
      )
      fitted_reduced <- rep(mean(means), length(means))
      residual_reduced <- means - fitted_reduced
      direct_statistics <- apply(plan$orders, 1L, function(index) {
        permuted <- fitted_reduced + residual_reduced[index]
        summary(aov(permuted ~ factor(unit_treatment)))[[1]][1L, "F value"]
      })
      extreme <- sum(direct_statistics >= direct_f - 1e-10)
      direct_p <- if (plan$exact) extreme / length(direct_statistics) else
        (extreme + 1) / (length(direct_statistics) + 1)
      cursor <- cursor + 1L
      rows[[cursor]] <- data.frame(
        scenario = scenario, iteration = iteration,
        engine_f = fit$test$statistic, engine_p = fit$test$p_value,
        direct_f = direct_f, direct_p = direct_p,
        max_permutation_difference = max(abs(
          fit$test$permutation_statistics - direct_statistics
        )), stringsAsFactors = FALSE
      )
    }
  }
  raw <- do.call(rbind, rows)
  summary <- do.call(rbind, lapply(names(scenarios), function(scenario) {
    selected <- raw$scenario == scenario
    do.call(rbind, lapply(c("engine", "direct"), function(method) {
      rejected <- sum(raw[[paste0(method, "_p")]][selected] <= alpha)
      uncertainty <- .binomial_summary(rejected, simulations)
      data.frame(
        scenario = scenario, method = method,
        estimand = if (scenario == "exchangeable_null") "type_I_error" else "power",
        simulations = simulations, rejected = rejected,
        rejection_rate = uncertainty[["rejection_rate"]],
        monte_carlo_se = uncertainty[["monte_carlo_se"]],
        confidence_low = uncertainty[["confidence_low"]],
        confidence_high = uncertainty[["confidence_high"]],
        stringsAsFactors = FALSE
      )
    }))
  }))
  rownames(summary) <- NULL
  out <- list(
    summary = summary,
    max_engine_direct_f_difference = max(abs(raw$engine_f - raw$direct_f)),
    max_engine_direct_p_difference = max(abs(raw$engine_p - raw$direct_p)),
    max_engine_direct_permutation_statistic_difference =
      max(raw$max_permutation_difference),
    results = if (keep_results) raw else NULL,
    settings = list(
      simulations = simulations, permutations = permutations,
      unit_sizes = sizes, effect = effect, cluster_sd = cluster_sd,
      residual_sd = residual_sd, alpha = alpha, seed = seed
    ),
    elapsed_seconds = proc.time()[["elapsed"]] - start
  )
  class(out) <- "permanova_unequal_random_intercept_validation"
  out
}

#' @export
print.permanova_unequal_random_intercept_validation <- function(x, ...) {
  cat("Unequal-replication random-intercept validation\n")
  print(x$summary, row.names = FALSE)
  cat("Maximum engine/direct F difference:",
      format(x$max_engine_direct_f_difference), "\n")
  cat("Maximum engine/direct permutation-statistic difference:",
      format(x$max_engine_direct_permutation_statistic_difference), "\n")
  invisible(x)
}

.phase12_nested_data <- function() {
  topology <- list(
    c(2L, 3L), c(2L, 3L), c(2L, 2L, 2L),
    c(2L, 2L, 2L), c(3L, 4L), c(3L, 4L),
    c(2L, 3L), c(2L, 3L), c(2L, 2L, 2L),
    c(2L, 2L, 2L), c(3L, 4L), c(3L, 4L)
  )
  treatment <- rep(c("control", "impact"), 6L)
  rows <- list()
  child_id <- 0L
  for (parent in seq_along(topology)) {
    for (size in topology[[parent]]) {
      child_id <- child_id + 1L
      rows[[length(rows) + 1L]] <- data.frame(
        region = parent, site = child_id, treatment = treatment[[parent]],
        replicate = seq_len(size)
      )
    }
  }
  data <- do.call(rbind, rows)
  data$region <- factor(data$region)
  data$site <- factor(data$site)
  data$treatment <- factor(data$treatment, levels = c("control", "impact"))
  data
}

.phase12_block_data <- function() {
  cells <- data.frame(
    block = factor(c(
      1, 1, 1, 2, 2, 2, 3, 3, 4, 4,
      5, 5, 5, 6, 6, 6, 7, 7, 8, 8
    )),
    treatment = factor(c(
      "A", "B", "C", "A", "B", "C", "A", "B", "B", "C",
      "A", "B", "C", "A", "B", "C", "A", "C", "A", "C"
    ), levels = c("A", "B", "C")),
    size = c(
      2L, 2L, 2L, 3L, 3L, 3L, 2L, 2L, 2L, 2L,
      4L, 4L, 4L, 2L, 2L, 2L, 3L, 3L, 3L, 3L
    )
  )
  data <- cells[rep(seq_len(nrow(cells)), cells$size),
                c("block", "treatment")]
  rownames(data) <- NULL
  data$replicate <- stats::ave(
    seq_len(nrow(data)), interaction(data$block, data$treatment),
    FUN = seq_along
  )
  data
}

#' Validate unbalanced nested and randomized-block inference
#'
#' Simulates the two Phase 12 envelopes under a null and location alternative,
#' comparing every distance statistic and p-value with an independent scalar
#' top-unit or cell-profile reconstruction.
#'
#' @param simulations Number of simulations per structure and scenario.
#' @param permutations Monte Carlo permutations per fit.
#' @param effect Treatment shift in the power scenario.
#' @param parent_sd,child_sd,block_sd,interaction_sd,residual_sd Generating
#'   standard deviations.
#' @param alpha Rejection threshold.
#' @param seed Reproducibility seed.
#' @param keep_results Retain iteration-level results.
#' @return A Phase 12 validation report.
#' @export
simulate_unbalanced_mixed_validation <- function(
    simulations = 500L, permutations = 99L, effect = 1,
    parent_sd = 1, child_sd = 0.7, block_sd = 1,
    interaction_sd = 1, residual_sd = 0.5, alpha = 0.05,
    seed = 20260719L, keep_results = FALSE) {
  simulations <- .positive_integer(simulations, "simulations")
  permutations <- .positive_integer(permutations, "permutations")
  scales <- c(parent_sd, child_sd, block_sd, interaction_sd, residual_sd)
  if (!is.finite(effect) || any(!is.finite(scales)) ||
      any(scales < 0) || residual_sd <= 0) {
    stop("Simulation effects must be finite and standard deviations valid.",
         call. = FALSE)
  }
  if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must lie strictly between zero and one.", call. = FALSE)
  }
  nested_data <- .phase12_nested_data()
  block_data <- .phase12_block_data()
  scenarios <- c(exchangeable_null = 0, location_shift = effect)
  rows <- vector("list", 2L * length(scenarios) * simulations)
  cursor <- 0L
  start <- proc.time()[["elapsed"]]

  for (scenario_index in seq_along(scenarios)) {
    scenario <- names(scenarios)[[scenario_index]]
    shift <- scenarios[[scenario_index]]
    for (iteration in seq_len(simulations)) {
      iteration_seed <- as.integer(seed + scenario_index * 100000L + iteration)

      set.seed(iteration_seed)
      nested_response <- shift * (nested_data$treatment == "impact") +
        rnorm(nlevels(nested_data$region), sd = parent_sd)[nested_data$region] +
        rnorm(nlevels(nested_data$site), sd = child_sd)[nested_data$site] +
        rnorm(nrow(nested_data), sd = residual_sd)
      nested_design <- compile_permanova_design(
        ~ treatment + (1 | region/site), nested_data, dist(nested_response)
      )
      nested_fit <- permanova_unbalanced_nested(
        nested_design, permutations = permutations,
        seed = iteration_seed + 5000L, keep_permutations = TRUE
      )
      nested_means <- aggregate(
        nested_response,
        list(region = nested_data$region,
             treatment = nested_data$treatment), mean
      )
      nested_means <- nested_means[
        order(as.integer(nested_means$region)), , drop = FALSE
      ]
      names(nested_means)[[3L]] <- "response"
      nested_direct_f <- summary(aov(
        response ~ treatment, nested_means
      ))[[1]][1L, "F value"]
      nested_structure <- .validate_unbalanced_nested(nested_design)
      set.seed(iteration_seed + 5000L)
      nested_plan <- .stratified_permutation_plan(
        nested_structure$topology_strata, permutations, 100000L
      )
      nested_reduced <- stats::lm(response ~ 1, nested_means)
      nested_statistics <- apply(nested_plan$orders, 1L, function(index) {
        permutation_data <- nested_means
        permutation_data$response <- stats::fitted(nested_reduced) +
          stats::residuals(nested_reduced)[index]
        summary(aov(response ~ treatment, permutation_data))[[1]][1L, "F value"]
      })
      nested_extreme <- sum(nested_statistics >= nested_direct_f - 1e-10)
      nested_direct_p <- if (nested_plan$exact) {
        nested_extreme / length(nested_statistics)
      } else (nested_extreme + 1) / (length(nested_statistics) + 1)
      cursor <- cursor + 1L
      rows[[cursor]] <- data.frame(
        structure = "nested", scenario = scenario, iteration = iteration,
        engine_f = nested_fit$test$statistic,
        engine_p = nested_fit$test$p_value,
        direct_f = nested_direct_f, direct_p = nested_direct_p,
        max_permutation_difference = max(abs(
          nested_fit$test$permutation_statistics - nested_statistics
        )), stringsAsFactors = FALSE
      )

      set.seed(iteration_seed + 200000L)
      block_response <- shift * (block_data$treatment != "A") +
        rnorm(nlevels(block_data$block), sd = block_sd)[block_data$block] +
        rnorm(nlevels(interaction(
          block_data$treatment, block_data$block, drop = TRUE
        )), sd = interaction_sd)[interaction(
          block_data$treatment, block_data$block, drop = TRUE
        )] + rnorm(nrow(block_data), sd = residual_sd)
      block_design <- compile_permanova_design(
        ~ treatment + (1 | block) + (1 | treatment:block),
        block_data, dist(block_response)
      )
      block_fit <- permanova_unbalanced_block(
        block_design, permutations = permutations,
        seed = iteration_seed + 205000L, keep_permutations = TRUE
      )
      block_structure <- .validate_unbalanced_block(block_design)
      block_means <- aggregate(
        block_response,
        list(treatment = block_data$treatment, block = block_data$block), mean
      )
      block_means <- block_means[match(
        paste(block_structure$cell_mapping$treatment,
              block_structure$cell_mapping$block),
        paste(block_means$treatment, block_means$block)
      ), , drop = FALSE]
      names(block_means)[[3L]] <- "response"
      block_table <- summary(aov(
        response ~ block + treatment, block_means
      ))[[1]]
      block_direct_f <- block_table["treatment", "Mean Sq"] /
        block_table["Residuals", "Mean Sq"]
      set.seed(iteration_seed + 205000L)
      block_plan <- .stratified_permutation_plan(
        block_structure$exchangeability_strata, permutations, 100000L
      )
      block_statistics <- apply(block_plan$orders, 1L, function(index) {
        permutation_data <- block_means
        permutation_data$response <- block_means$response[index]
        table <- summary(aov(
          response ~ block + treatment, permutation_data
        ))[[1]]
        table["treatment", "Mean Sq"] / table["Residuals", "Mean Sq"]
      })
      block_extreme <- sum(block_statistics >= block_direct_f - 1e-10)
      block_direct_p <- if (block_plan$exact) {
        block_extreme / length(block_statistics)
      } else (block_extreme + 1) / (length(block_statistics) + 1)
      cursor <- cursor + 1L
      rows[[cursor]] <- data.frame(
        structure = "block", scenario = scenario, iteration = iteration,
        engine_f = block_fit$test$statistic,
        engine_p = block_fit$test$p_value,
        direct_f = block_direct_f, direct_p = block_direct_p,
        max_permutation_difference = max(abs(
          block_fit$test$permutation_statistics - block_statistics
        )), stringsAsFactors = FALSE
      )
    }
  }
  raw <- do.call(rbind, rows)
  summary <- do.call(rbind, lapply(c("nested", "block"), function(structure) {
    do.call(rbind, lapply(names(scenarios), function(scenario) {
      selected <- raw$structure == structure & raw$scenario == scenario
      do.call(rbind, lapply(c("engine", "direct"), function(method) {
        rejected <- sum(raw[[paste0(method, "_p")]][selected] <= alpha)
        uncertainty <- .binomial_summary(rejected, simulations)
        data.frame(
          structure = structure, scenario = scenario, method = method,
          estimand = if (scenario == "exchangeable_null") "type_I_error" else "power",
          simulations = simulations, rejected = rejected,
          rejection_rate = uncertainty[["rejection_rate"]],
          monte_carlo_se = uncertainty[["monte_carlo_se"]],
          confidence_low = uncertainty[["confidence_low"]],
          confidence_high = uncertainty[["confidence_high"]],
          stringsAsFactors = FALSE
        )
      }))
    }))
  }))
  rownames(summary) <- NULL
  out <- list(
    summary = summary,
    max_engine_direct_f_difference = max(abs(raw$engine_f - raw$direct_f)),
    max_engine_direct_p_difference = max(abs(raw$engine_p - raw$direct_p)),
    max_engine_direct_permutation_statistic_difference =
      max(raw$max_permutation_difference),
    results = if (keep_results) raw else NULL,
    settings = list(
      simulations = simulations, permutations = permutations, effect = effect,
      parent_sd = parent_sd, child_sd = child_sd, block_sd = block_sd,
      interaction_sd = interaction_sd, residual_sd = residual_sd,
      alpha = alpha, seed = seed
    ),
    elapsed_seconds = proc.time()[["elapsed"]] - start
  )
  class(out) <- "permanova_unbalanced_mixed_validation"
  out
}

#' @export
print.permanova_unbalanced_mixed_validation <- function(x, ...) {
  cat("Unbalanced nested and randomized-block validation\n")
  print(x$summary, row.names = FALSE)
  cat("Maximum engine/direct F difference:",
      format(x$max_engine_direct_f_difference), "\n")
  cat("Maximum engine/direct permutation-statistic difference:",
      format(x$max_engine_direct_permutation_statistic_difference), "\n")
  invisible(x)
}

#' Validate planned one-factor contrasts
#'
#' Simulates unequal group counts and compares a planned A-minus-B distance
#' contrast with an independently constrained scalar model and the classical
#' one-df contrast F test under a null and location alternative.
#'
#' @param simulations Number of simulations per scenario.
#' @param permutations Monte Carlo permutations per fit.
#' @param effect B-versus-A location shift in the power scenario.
#' @param residual_sd Generating residual standard deviation.
#' @param alpha Rejection threshold.
#' @param seed Reproducibility seed.
#' @param keep_results Retain iteration-level results.
#' @return A Phase 13 validation report.
#' @export
simulate_contrast_validation <- function(
    simulations = 500L, permutations = 99L, effect = 0.8,
    residual_sd = 1, alpha = 0.05, seed = 20260719L,
    keep_results = FALSE) {
  simulations <- .positive_integer(simulations, "simulations")
  permutations <- .positive_integer(permutations, "permutations")
  if (!is.finite(effect) || !is.finite(residual_sd) || residual_sd <= 0) {
    stop("`effect` must be finite and `residual_sd` positive.", call. = FALSE)
  }
  if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must lie strictly between zero and one.", call. = FALSE)
  }
  counts <- c(A = 4L, B = 6L, C = 5L)
  data <- data.frame(
    group = factor(rep(names(counts), counts), levels = names(counts))
  )
  coefficients <- c(A = 1, B = -1, C = 0)
  X_group <- stats::model.matrix(~ 0 + group, data)
  Q <- qr.Q(qr(matrix(coefficients, ncol = 1L)), complete = TRUE)
  X_reduced <- X_group %*% Q[, -1L, drop = FALSE]
  scenarios <- c(exchangeable_null = 0, location_shift = effect)
  rows <- vector("list", length(scenarios) * simulations)
  cursor <- 0L
  start <- proc.time()[["elapsed"]]
  for (scenario_index in seq_along(scenarios)) {
    scenario <- names(scenarios)[[scenario_index]]
    shift <- scenarios[[scenario_index]]
    for (iteration in seq_len(simulations)) {
      iteration_seed <- as.integer(seed + scenario_index * 100000L + iteration)
      set.seed(iteration_seed)
      response <- shift * (data$group == "B") +
        stats::rnorm(nrow(data), sd = residual_sd)
      design <- compile_permanova_design(~ group, data, dist(response))
      fit <- permanova_contrasts(
        design, coefficients, adjust = "none",
        permutations = permutations, seed = iteration_seed + 5000L,
        keep_permutations = TRUE
      )
      group_means <- tapply(response, data$group, mean)
      contrast_ss <- sum(coefficients * group_means) ^ 2 /
        sum(coefficients ^ 2 / counts)
      full <- stats::lm(response ~ group, data)
      full_ms <- sum(stats::residuals(full) ^ 2) / stats::df.residual(full)
      direct_f <- contrast_ss / full_ms
      classical_p <- stats::pf(
        direct_f, 1L, stats::df.residual(full), lower.tail = FALSE
      )

      reduced <- stats::lm.fit(X_reduced, response)
      fitted_reduced <- as.vector(X_reduced %*% reduced$coefficients)
      residual_reduced <- response - fitted_reduced
      set.seed(iteration_seed + 5000L)
      plan <- .permutation_orders(
        nrow(data), permutations = permutations, max_exact = 100000L
      )
      direct_statistics <- apply(plan$orders, 1L, function(index) {
        permuted <- fitted_reduced + residual_reduced[index]
        means <- tapply(permuted, data$group, mean)
        ss <- sum(coefficients * means) ^ 2 /
          sum(coefficients ^ 2 / counts)
        permuted_full <- stats::lm(permuted ~ group, data)
        ss / (sum(stats::residuals(permuted_full) ^ 2) /
                stats::df.residual(permuted_full))
      })
      extreme <- sum(direct_statistics >= direct_f - 1e-10)
      direct_p <- if (plan$exact) extreme / length(direct_statistics) else
        (extreme + 1) / (length(direct_statistics) + 1)
      cursor <- cursor + 1L
      rows[[cursor]] <- data.frame(
        scenario = scenario, iteration = iteration,
        engine_f = fit$tests[[1L]]$statistic,
        engine_p = fit$tests[[1L]]$p_value,
        direct_f = direct_f, direct_p = direct_p,
        classical_f = direct_f, classical_p = classical_p,
        max_permutation_difference = max(abs(
          fit$tests[[1L]]$permutation_statistics - direct_statistics
        )), stringsAsFactors = FALSE
      )
    }
  }
  raw <- do.call(rbind, rows)
  summary <- do.call(rbind, lapply(names(scenarios), function(scenario) {
    selected <- raw$scenario == scenario
    do.call(rbind, lapply(c("engine", "direct", "classical"), function(method) {
      rejected <- sum(raw[[paste0(method, "_p")]][selected] <= alpha)
      uncertainty <- .binomial_summary(rejected, simulations)
      data.frame(
        scenario = scenario, method = method,
        estimand = if (scenario == "exchangeable_null") "type_I_error" else "power",
        simulations = simulations, rejected = rejected,
        rejection_rate = uncertainty[["rejection_rate"]],
        monte_carlo_se = uncertainty[["monte_carlo_se"]],
        confidence_low = uncertainty[["confidence_low"]],
        confidence_high = uncertainty[["confidence_high"]],
        stringsAsFactors = FALSE
      )
    }))
  }))
  rownames(summary) <- NULL
  out <- list(
    summary = summary,
    max_engine_direct_f_difference = max(abs(raw$engine_f - raw$direct_f)),
    max_engine_direct_p_difference = max(abs(raw$engine_p - raw$direct_p)),
    max_engine_direct_permutation_statistic_difference =
      max(raw$max_permutation_difference),
    max_engine_classical_f_difference =
      max(abs(raw$engine_f - raw$classical_f)),
    results = if (keep_results) raw else NULL,
    settings = list(
      simulations = simulations, permutations = permutations,
      group_counts = counts, coefficients = coefficients,
      effect = effect, residual_sd = residual_sd,
      alpha = alpha, seed = seed
    ),
    elapsed_seconds = proc.time()[["elapsed"]] - start
  )
  class(out) <- "permanova_contrast_validation"
  out
}

#' @export
print.permanova_contrast_validation <- function(x, ...) {
  cat("Planned distance-contrast validation\n")
  print(x$summary, row.names = FALSE)
  cat("Maximum engine/direct F difference:",
      format(x$max_engine_direct_f_difference), "\n")
  cat("Maximum engine/direct permutation-statistic difference:",
      format(x$max_engine_direct_permutation_statistic_difference), "\n")
  invisible(x)
}

.audit_assignment_matrix <- function(impact_per_stratum) {
  impact_per_stratum <- as.integer(impact_per_stratum)
  if (length(impact_per_stratum) != 3L ||
      any(impact_per_stratum < 1L | impact_per_stratum > 3L)) {
    stop("Audit allocations require 1--3 impact units in each four-unit size stratum.",
         call. = FALSE)
  }
  sets <- lapply(impact_per_stratum, function(k) {
    utils::combn(4L, k, simplify = FALSE)
  })
  grid <- expand.grid(lapply(sets, seq_along))
  assignments <- matrix(FALSE, nrow(grid), 12L)
  for (row in seq_len(nrow(grid))) {
    for (stratum in seq_len(3L)) {
      positions <- (stratum - 1L) * 4L +
        sets[[stratum]][[grid[row, stratum]]]
      assignments[row, positions] <- TRUE
    }
  }
  assignments
}

.audit_two_group_f <- function(centroids, assignments) {
  impact_n <- rowSums(assignments)
  control_n <- length(centroids) - impact_n
  impact_mean <- as.vector(assignments %*% centroids) / impact_n
  control_mean <- as.vector((!assignments) %*% centroids) / control_n
  grand_mean <- mean(centroids)
  between <- impact_n * (impact_mean - grand_mean) ^ 2 +
    control_n * (control_mean - grand_mean) ^ 2
  total <- sum((centroids - grand_mean) ^ 2)
  between / ((total - between) / (length(centroids) - 2L))
}

.audit_null_response <- function(distribution, data, unit_treatment,
                                 cluster_sd, residual_sd) {
  if (distribution == "heavy_tail") {
    unit_effect <- stats::rt(nlevels(data$site), df = 3) /
      sqrt(3) * cluster_sd
    residual <- stats::rt(nrow(data), df = 3) /
      sqrt(3) * residual_sd
  } else {
    unit_effect <- stats::rnorm(nlevels(data$site), sd = cluster_sd)
    residual_scale <- residual_sd * if (distribution == "dispersion_violation") {
      ifelse(unit_treatment[data$site] == "impact", 2, 1)
    } else 1
    residual <- stats::rnorm(nrow(data), sd = residual_scale)
  }
  unit_effect[data$site] + residual
}

#' Audit unequal-random-intercept calibration across independent batches
#'
#' Uses exact unique treatment assignments within cluster-size strata and a
#' separate Monte Carlo approximation. Batch seeds are separated by one million
#' to prevent the overlapping per-iteration streams that make adjacent-seed
#' summaries appear artificially stable. Gaussian, heavy-tail, unequal-
#' allocation, and deliberate dispersion-violation scenarios are included.
#'
#' @param batches Number of independent seed batches.
#' @param simulations_per_batch Null simulations in each batch and scenario.
#' @param permutations Monte Carlo assignment draws per simulation.
#' @param engine_checks Number of first-batch Gaussian datasets also passed
#'   through the complete distance engine with exhaustive permutations.
#' @param cluster_sd,residual_sd Generating standard deviations.
#' @param alpha Rejection threshold.
#' @param seed Base reproducibility seed. Batch seeds are spaced internally.
#' @param keep_results Retain simulation-level p-values.
#' @return A Phase 14 exact-randomization calibration audit.
#' @export
audit_unequal_random_intercept_calibration <- function(
    batches = 20L, simulations_per_batch = 500L,
    permutations = 99L, engine_checks = 10L,
    cluster_sd = 1, residual_sd = 1, alpha = 0.05,
    seed = 20260719L, keep_results = FALSE) {
  batches <- .positive_integer(batches, "batches")
  simulations_per_batch <- .positive_integer(
    simulations_per_batch, "simulations_per_batch"
  )
  permutations <- .positive_integer(permutations, "permutations")
  if (length(engine_checks) != 1L || !is.finite(engine_checks) ||
      engine_checks < 0 || engine_checks != as.integer(engine_checks)) {
    stop("`engine_checks` must be one non-negative integer.", call. = FALSE)
  }
  engine_checks <- as.integer(engine_checks)
  if (engine_checks > simulations_per_batch) {
    stop("`engine_checks` cannot exceed `simulations_per_batch`.",
         call. = FALSE)
  }
  if (simulations_per_batch >= 50000L) {
    stop("Use fewer than 50000 simulations per batch to retain disjoint seed ranges.",
         call. = FALSE)
  }
  if (any(!is.finite(c(cluster_sd, residual_sd))) || cluster_sd < 0 ||
      residual_sd <= 0 || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("Simulation scales must be valid and `alpha` strictly between zero and one.",
         call. = FALSE)
  }
  batch_seeds <- as.double(seed) + (seq_len(batches) - 1L) * 1000000
  if (any(batch_seeds > .Machine$integer.max - 500000L) ||
      any(batch_seeds < 0)) {
    stop("The requested base seed and batch count exceed supported integer seed ranges.",
         call. = FALSE)
  }
  batch_seeds <- as.integer(batch_seeds)
  scenarios <- data.frame(
    scenario = c(
      "gaussian_balanced", "heavy_tail_balanced",
      "gaussian_unequal_allocation", "dispersion_violation"
    ),
    distribution = c("gaussian", "heavy_tail", "gaussian",
                     "dispersion_violation"),
    allocation = c("2,2,2", "2,2,2", "1,2,3", "2,2,2"),
    exchangeability = c(TRUE, TRUE, TRUE, FALSE),
    stringsAsFactors = FALSE
  )
  sizes <- rep(c(2L, 3L, 5L), each = 4L)
  site <- factor(rep(seq_along(sizes), sizes))
  data_template <- data.frame(site = site)
  rows <- vector("list", nrow(scenarios) * batches * simulations_per_batch)
  cursor <- 0L
  max_engine_f_difference <- 0
  max_engine_exact_p_difference <- 0
  start <- proc.time()[["elapsed"]]

  for (scenario_index in seq_len(nrow(scenarios))) {
    allocation <- as.integer(strsplit(
      scenarios$allocation[[scenario_index]], ",", fixed = TRUE
    )[[1L]])
    assignments <- .audit_assignment_matrix(allocation)
    observed <- logical(12L)
    for (stratum in seq_len(3L)) {
      observed[(stratum - 1L) * 4L + seq_len(allocation[[stratum]])] <- TRUE
    }
    observed_row <- which(apply(assignments, 1L, function(x) {
      all(x == observed)
    }))
    unit_treatment <- factor(
      ifelse(observed, "impact", "control"),
      levels = c("control", "impact")
    )
    engine_data <- data.frame(
      site = site,
      treatment = factor(unit_treatment[site],
                         levels = c("control", "impact"))
    )

    for (batch in seq_len(batches)) {
      for (iteration in seq_len(simulations_per_batch)) {
        iteration_seed <- as.integer(
          batch_seeds[[batch]] + scenario_index * 100000L + iteration
        )
        set.seed(iteration_seed)
        response <- .audit_null_response(
          scenarios$distribution[[scenario_index]], data_template,
          unit_treatment, cluster_sd, residual_sd
        )
        centroids <- as.numeric(tapply(response, site, mean))
        statistics <- .audit_two_group_f(centroids, assignments)
        observed_f <- statistics[[observed_row]]
        exact_p <- mean(statistics >= observed_f - 1e-10)
        set.seed(iteration_seed + 5000L)
        sampled <- sample.int(nrow(assignments), permutations, replace = TRUE)
        monte_carlo_p <- (
          sum(statistics[sampled] >= observed_f - 1e-10) + 1
        ) / (permutations + 1)

        if (scenario_index == 1L && batch == 1L &&
            iteration <= engine_checks) {
          engine_design <- compile_permanova_design(
            ~ treatment + (1 | site), engine_data, dist(response)
          )
          engine <- suppressWarnings(permanova_unequal_random_intercept(
            engine_design, permutations = "all", max_exact = 20000L
          ))
          max_engine_f_difference <- max(
            max_engine_f_difference,
            abs(engine$test$statistic - observed_f)
          )
          max_engine_exact_p_difference <- max(
            max_engine_exact_p_difference,
            abs(engine$test$p_value - exact_p)
          )
        }
        cursor <- cursor + 1L
        rows[[cursor]] <- data.frame(
          scenario = scenarios$scenario[[scenario_index]],
          exchangeability = scenarios$exchangeability[[scenario_index]],
          batch = batch, batch_seed = batch_seeds[[batch]],
          iteration = iteration, observed_f = observed_f,
          exact_p = exact_p, monte_carlo_p = monte_carlo_p,
          exact_assignments = nrow(assignments),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  raw <- do.call(rbind, rows)
  batch_summary <- do.call(rbind, lapply(
    split(raw, list(raw$scenario, raw$batch), drop = TRUE), function(part) {
      data.frame(
        scenario = part$scenario[[1L]],
        exchangeability = part$exchangeability[[1L]],
        batch = part$batch[[1L]], batch_seed = part$batch_seed[[1L]],
        simulations = nrow(part),
        exact_rejection_rate = mean(part$exact_p <= alpha),
        monte_carlo_rejection_rate = mean(part$monte_carlo_p <= alpha),
        stringsAsFactors = FALSE
      )
    }
  ))
  rownames(batch_summary) <- NULL
  summary <- do.call(rbind, lapply(
    split(raw, raw$scenario), function(part) {
      do.call(rbind, lapply(c("exact", "monte_carlo"), function(method) {
        p_values <- part[[paste0(method, "_p")]]
        rejected <- sum(p_values <= alpha)
        uncertainty <- .binomial_summary(rejected, length(p_values))
        data.frame(
          scenario = part$scenario[[1L]],
          exchangeability = part$exchangeability[[1L]], method = method,
          simulations = length(p_values), rejected = rejected,
          rejection_rate = uncertainty[["rejection_rate"]],
          monte_carlo_se = uncertainty[["monte_carlo_se"]],
          confidence_low = uncertainty[["confidence_low"]],
          confidence_high = uncertainty[["confidence_high"]],
          stringsAsFactors = FALSE
        )
      }))
    }
  ))
  rownames(summary) <- NULL
  valid_exact <- summary$exchangeability & summary$method == "exact"
  calibration_pass <- all(summary$confidence_low[valid_exact] <= alpha)
  out <- list(
    summary = summary, batch_summary = batch_summary,
    calibration_pass = calibration_pass,
    pinned_batch = batch_summary[batch_summary$batch == 1L, , drop = FALSE],
    mean_absolute_monte_carlo_exact_p_difference =
      mean(abs(raw$monte_carlo_p - raw$exact_p)),
    max_engine_exact_f_difference = max_engine_f_difference,
    max_engine_exact_p_difference = max_engine_exact_p_difference,
    seed_design = list(
      base_seed = seed, batch_seeds = batch_seeds,
      spacing = 1000000L,
      note = paste(
        "Non-overlapping batch seeds are required because each batch assigns",
        "deterministic per-iteration seeds."
      )
    ),
    settings = list(
      batches = batches, simulations_per_batch = simulations_per_batch,
      permutations = permutations, engine_checks = engine_checks,
      cluster_sd = cluster_sd, residual_sd = residual_sd,
      alpha = alpha, seed = seed
    ),
    results = if (keep_results) raw else NULL,
    elapsed_seconds = proc.time()[["elapsed"]] - start
  )
  class(out) <- "permanova_random_intercept_calibration_audit"
  out
}

#' @export
print.permanova_random_intercept_calibration_audit <- function(x, ...) {
  cat("Unequal-random-intercept exact calibration audit\n")
  print(x$summary, row.names = FALSE)
  cat("Valid-scenario calibration gate:",
      if (x$calibration_pass) "PASS" else "FAIL", "\n")
  cat("Maximum full-engine/exact F difference:",
      format(x$max_engine_exact_f_difference), "\n")
  cat("Maximum full-engine/exact p-value difference:",
      format(x$max_engine_exact_p_difference), "\n")
  invisible(x)
}

#' Validate multifactor marginal contrasts
#'
#' Simulates an unequal-cell two-factor design and compares an equal-weight
#' marginal A contrast with a direct scalar general linear hypothesis and
#' response-level Freedman--Lane reconstruction.
#'
#' @param simulations Number of simulations per scenario.
#' @param permutations Monte Carlo permutations per fit.
#' @param effect Equal-grid marginal A shift in the power scenario.
#' @param residual_sd Generating residual standard deviation.
#' @param alpha Rejection threshold.
#' @param seed Reproducibility seed.
#' @param keep_results Retain iteration-level results.
#' @return A Phase 15 validation report.
#' @export
simulate_marginal_contrast_validation <- function(
    simulations = 500L, permutations = 99L, effect = 0.8,
    residual_sd = 1, alpha = 0.05, seed = 20260719L,
    keep_results = FALSE) {
  simulations <- .positive_integer(simulations, "simulations")
  permutations <- .positive_integer(permutations, "permutations")
  if (!is.finite(effect) || !is.finite(residual_sd) || residual_sd <= 0 ||
      !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("Simulation settings must be finite with positive scale and valid alpha.",
         call. = FALSE)
  }
  cells <- expand.grid(
    A = factor(c("a1", "a2")), B = factor(c("b1", "b2", "b3"))
  )
  counts <- c(3L, 5L, 4L, 2L, 6L, 3L)
  data <- cells[rep(seq_len(nrow(cells)), counts), , drop = FALSE]
  rownames(data) <- NULL
  B_effect <- c(b1 = -0.6, b2 = 0.1, b3 = 0.7)
  interaction <- c(b1 = -0.5, b2 = 0, b3 = 0.5)
  contrast <- c(a1 = -1, a2 = 1)
  scenarios <- c(exchangeable_null = 0, location_shift = effect)
  rows <- vector("list", length(scenarios) * simulations)
  cursor <- 0L
  start <- proc.time()[["elapsed"]]
  for (scenario_index in seq_along(scenarios)) {
    scenario <- names(scenarios)[[scenario_index]]
    shift <- scenarios[[scenario_index]]
    for (iteration in seq_len(simulations)) {
      iteration_seed <- as.integer(seed + scenario_index * 100000L + iteration)
      set.seed(iteration_seed)
      response <- B_effect[data$B] +
        (data$A == "a2") * (shift + interaction[data$B]) +
        stats::rnorm(nrow(data), sd = residual_sd)
      design <- compile_permanova_design(~ A * B, data, dist(response))
      fit <- permanova_marginal_contrasts(
        design, "A", contrast, weights = "equal", adjust = "none",
        permutations = permutations, seed = iteration_seed + 5000L,
        keep_permutations = TRUE
      )
      X <- design$X
      coefficient_row <- fit$tests[[1L]]$marginal_coefficient_row
      full <- stats::lm.fit(X, response)
      information_inverse <- solve(crossprod(X))
      estimate <- as.numeric(coefficient_row %*% full$coefficients)
      contrast_variance <- as.numeric(
        coefficient_row %*% information_inverse %*% coefficient_row
      )
      contrast_ss <- estimate ^ 2 / contrast_variance
      residual_ms <- sum(full$residuals ^ 2) / full$df.residual
      direct_f <- contrast_ss / residual_ms
      classical_p <- stats::pf(
        direct_f, 1L, full$df.residual, lower.tail = FALSE
      )
      projector <- .linear_contrast_projector(
        X, coefficient_row, projection_matrix(X)$H, 1e-10
      )
      fitted_reduced <- as.vector(projector$H_reduced %*% response)
      residual_reduced <- response - fitted_reduced
      set.seed(iteration_seed + 5000L)
      plan <- .permutation_orders(nrow(data), permutations, 100000L)
      direct_statistics <- apply(plan$orders, 1L, function(index) {
        y <- fitted_reduced + residual_reduced[index]
        model <- stats::lm.fit(X, y)
        estimate <- as.numeric(coefficient_row %*% model$coefficients)
        ss <- estimate ^ 2 / contrast_variance
        ss / (sum(model$residuals ^ 2) / model$df.residual)
      })
      extreme <- sum(direct_statistics >= direct_f - 1e-10)
      direct_p <- if (plan$exact) extreme / length(direct_statistics) else
        (extreme + 1) / (length(direct_statistics) + 1)
      cursor <- cursor + 1L
      rows[[cursor]] <- data.frame(
        scenario = scenario, iteration = iteration,
        engine_f = fit$tests[[1L]]$statistic,
        engine_p = fit$tests[[1L]]$p_value,
        direct_f = direct_f, direct_p = direct_p,
        classical_f = direct_f, classical_p = classical_p,
        max_permutation_difference = max(abs(
          fit$tests[[1L]]$permutation_statistics - direct_statistics
        )), stringsAsFactors = FALSE
      )
    }
  }
  raw <- do.call(rbind, rows)
  summary <- do.call(rbind, lapply(names(scenarios), function(scenario) {
    selected <- raw$scenario == scenario
    do.call(rbind, lapply(c("engine", "direct", "classical"), function(method) {
      rejected <- sum(raw[[paste0(method, "_p")]][selected] <= alpha)
      uncertainty <- .binomial_summary(rejected, simulations)
      data.frame(
        scenario = scenario, method = method,
        estimand = if (scenario == "exchangeable_null") "type_I_error" else "power",
        simulations = simulations, rejected = rejected,
        rejection_rate = uncertainty[["rejection_rate"]],
        monte_carlo_se = uncertainty[["monte_carlo_se"]],
        confidence_low = uncertainty[["confidence_low"]],
        confidence_high = uncertainty[["confidence_high"]],
        stringsAsFactors = FALSE
      )
    }))
  }))
  rownames(summary) <- NULL
  out <- list(
    summary = summary,
    max_engine_direct_f_difference = max(abs(raw$engine_f - raw$direct_f)),
    max_engine_direct_p_difference = max(abs(raw$engine_p - raw$direct_p)),
    max_engine_direct_permutation_statistic_difference =
      max(raw$max_permutation_difference),
    max_engine_classical_f_difference =
      max(abs(raw$engine_f - raw$classical_f)),
    results = if (keep_results) raw else NULL,
    settings = list(
      simulations = simulations, permutations = permutations,
      cell_counts = counts, weighting = "equal", effect = effect,
      residual_sd = residual_sd, alpha = alpha, seed = seed
    ),
    elapsed_seconds = proc.time()[["elapsed"]] - start
  )
  class(out) <- "permanova_marginal_contrast_validation"
  out
}

#' @export
print.permanova_marginal_contrast_validation <- function(x, ...) {
  cat("Multifactor marginal-contrast validation\n")
  print(x$summary, row.names = FALSE)
  cat("Maximum engine/direct F difference:",
      format(x$max_engine_direct_f_difference), "\n")
  cat("Maximum engine/direct permutation-statistic difference:",
      format(x$max_engine_direct_permutation_statistic_difference), "\n")
  invisible(x)
}

.validation_row <- function(check, passed, metric, tolerance, detail) {
  data.frame(check = check, passed = isTRUE(passed), metric = metric,
             tolerance = tolerance, detail = detail, stringsAsFactors = FALSE)
}

#' Run the deterministic validation suite
#'
#' @param seed Reproducibility seed.
#' @param include_vegan Compare the supported fixed one-way case with
#'   `vegan::adonis2()` when vegan is installed.
#' @return A data frame of validation checks with class `permanova_validation`.
#' @export
run_permanova_validation <- function(seed = 1L, include_vegan = TRUE) {
  set.seed(seed)
  rows <- list()
  y <- c(1.1, 1.5, 1.7, 4.2, 4.8, 5.1)
  group <- factor(rep(c("A", "B"), each = 3))
  one <- permanova_oneway(dist(y), group, permutations = 0, engine = "projection")
  reference <- unname(summary(aov(y ~ group))[[1]][["F value"]][1])
  difference <- abs(one$table$pseudo_f[[1]] - reference)
  rows[[length(rows) + 1L]] <- .validation_row(
    "one-way Euclidean F equals ANOVA", difference < 1e-10,
    difference, 1e-10, "analytical"
  )

  generated <- simulate_balanced_mixed(
    groups_per_treatment = 3, replicates = 3, effect = 1, seed = seed
  )
  mixed_design <- compile_permanova_design(
    ~ treatment + (1 | site), generated$data, generated$distance
  )
  mixed <- suppressWarnings(permanova_mixed(mixed_design, permutations = 0))
  means <- aggregate(generated$response[, 1L],
                     list(treatment = generated$data$treatment,
                          site = generated$data$site), mean)
  names(means)[[3L]] <- "response"
  mixed_reference <- unname(summary(aov(response ~ treatment, means))[[1]][["F value"]][1])
  mixed_difference <- abs(mixed$tests$treatment$statistic - mixed_reference)
  rows[[length(rows) + 1L]] <- .validation_row(
    "mixed Euclidean F equals whole-plot ANOVA", mixed_difference < 1e-10,
    mixed_difference, 1e-10, "analytical EMS denominator"
  )

  small <- simulate_balanced_mixed(groups_per_treatment = 2, replicates = 2, seed = seed)
  small_design <- compile_permanova_design(
    ~ treatment + (1 | site), small$data, small$distance
  )
  exact <- suppressWarnings(permanova_mixed(small_design, permutations = "all"))
  exact_ok <- exact$permutation$attainable == 24 && exact$permutation$used == 24 &&
    abs(exact$tests$treatment$p_value * 24 -
          round(exact$tests$treatment$p_value * 24)) < 1e-12
  rows[[length(rows) + 1L]] <- .validation_row(
    "exhaustive whole-plot enumeration", exact_ok,
    exact$permutation$used, 24, "all 4! arrangements"
  )

  semimetric <- matrix(c(0, 1, 1, 1, 0, 3, 1, 3, 0), 3, 3)
  spectrum <- distance_diagnostics(semimetric)
  rows[[length(rows) + 1L]] <- .validation_row(
    "negative eigenvalues retained", !spectrum$euclidean && spectrum$negative_count > 0,
    spectrum$negative_count, 1, "non-Euclidean diagnostic"
  )

  unbalanced_data <- generated$data[-1L, ]
  unbalanced_distance <- as.matrix(generated$distance)[-1L, -1L]
  unbalanced_design <- compile_permanova_design(
    ~ treatment + (1 | site), unbalanced_data, unbalanced_distance
  )
  rejection <- tryCatch({
    permanova_mixed(unbalanced_design, permutations = 0)
    "not rejected"
  }, error = conditionMessage)
  rows[[length(rows) + 1L]] <- .validation_row(
    "unsupported imbalance rejected", grepl("equal observations", rejection),
    as.numeric(grepl("equal observations", rejection)), 1, rejection
  )

  x <- seq(-1, 1, length.out = 20)
  projection <- projection_matrix(cbind(1, x, x + 1e-12))
  idempotence_error <- max(abs(projection$H %*% projection$H - projection$H))
  rows[[length(rows) + 1L]] <- .validation_row(
    "rank-safe projection stability", idempotence_error < 1e-10,
    idempotence_error, 1e-10, "near-collinear matrix"
  )

  crossed_data <- expand.grid(
    replicate = factor(1:2), treatment = factor(c("control", "impact")),
    block = factor(1:4)
  )
  block_effect <- c(-0.6, -0.1, 0.2, 0.5)
  interaction_effect <- c(-1.1, 0.4, 0.8, -0.2, 0.9, -0.7, 0.3, -0.5)
  crossed_cell <- interaction(crossed_data$treatment, crossed_data$block, drop = TRUE)
  crossed_response <- block_effect[crossed_data$block] +
    interaction_effect[crossed_cell] + rep(c(-0.1, 0.1), 8)
  crossed_design <- compile_permanova_design(
    ~ treatment + (1 | block) + (1 | treatment:block),
    crossed_data, dist(crossed_response)
  )
  crossed_fit <- permanova_crossed(crossed_design)
  crossed_reference <- summary(aov(
    crossed_response ~ treatment * block, crossed_data
  ))[[1]]
  crossed_rows <- trimws(rownames(crossed_reference))
  crossed_reference_f <- crossed_reference[["Mean Sq"]][match("treatment", crossed_rows)] /
    crossed_reference[["Mean Sq"]][match("treatment:block", crossed_rows)]
  crossed_difference <- abs(
    crossed_fit$tests$treatment$statistic - crossed_reference_f
  )
  rows[[length(rows) + 1L]] <- .validation_row(
    "declared fixed-by-random denominator equals classical interaction test",
    crossed_difference < 1e-10, crossed_difference, 1e-10,
    crossed_fit$tests$treatment$denominator$expression
  )

  crossed_permutation <- suppressWarnings(permanova_crossed(
    crossed_design, p_method = "permutation", permutations = "all",
    keep_permutations = TRUE
  ))
  crossed_permutation_ok <- crossed_permutation$permutation$exact &&
    crossed_permutation$permutation$attainable == 2 ^ 4 &&
    crossed_permutation$permutation$used == 2 ^ 4 &&
    crossed_permutation$tests$treatment$unique_permutations_used == 2 ^ 4
  rows[[length(rows) + 1L]] <- .validation_row(
    "exhaustive crossed within-block cell enumeration",
    crossed_permutation_ok, crossed_permutation$permutation$used, 2 ^ 4,
    "all 2^4 intact treatment-cell arrangements"
  )

  crossed_structure <- .validate_crossed_design(crossed_design)
  crossed_permutation_structure <- .validate_crossed_permutation_design(
    crossed_design, crossed_structure
  )
  crossed_plan <- .crossed_block_permutation_plan(
    crossed_permutation_structure, "all", 100L
  )
  crossed_direct_statistics <- vapply(
    seq_len(nrow(crossed_plan$orders)),
    function(row) {
      index <- .block_permutation_index(
        crossed_permutation_structure$cell_blocks,
        crossed_plan$orders[row, ]
      )
      reference <- summary(aov(
        crossed_response[index] ~ treatment * block, crossed_data
      ))[[1]]
      reference_rows <- trimws(rownames(reference))
      reference[["Mean Sq"]][match("treatment", reference_rows)] /
        reference[["Mean Sq"]][match("treatment:block", reference_rows)]
    },
    numeric(1)
  )
  crossed_permutation_difference <- max(abs(
    crossed_permutation$tests$treatment$permutation_statistics -
      crossed_direct_statistics
  ))
  rows[[length(rows) + 1L]] <- .validation_row(
    "crossed permutation statistics equal direct blocked ANOVA",
    crossed_permutation_difference < 1e-10,
    crossed_permutation_difference, 1e-10,
    "response-level reference for all arrangements"
  )

  provenance_ok <- nrow(crossed_design$containment_graph) == 0L &&
    crossed_design$random_term_table$origin[[2L]] == "interaction"
  rows[[length(rows) + 1L]] <- .validation_row(
    "explicit colon interaction retains non-nested provenance", provenance_ok,
    as.numeric(provenance_ok), 1, crossed_design$random_term_table$origin[[2L]]
  )

  additive_design <- compile_permanova_design(
    ~ treatment + (1 | block), crossed_data, dist(crossed_response)
  )
  additive_message <- tryCatch({
    permanova_crossed(additive_design)
    "not rejected"
  }, error = conditionMessage)
  additive_rejected <- grepl("Unmodelled fixed-by-random", additive_message)
  rows[[length(rows) + 1L]] <- .validation_row(
    "omitted fixed-by-random interaction rejected by default", additive_rejected,
    as.numeric(additive_rejected), 1, additive_message
  )

  non_euclidean <- as.matrix(dist(cbind(seq_len(nrow(crossed_data)), crossed_response)))
  non_euclidean[1, 2] <- non_euclidean[2, 1] <- 4 * max(non_euclidean)
  non_euclidean_design <- compile_permanova_design(
    ~ treatment + (1 | block) + (1 | treatment:block),
    crossed_data, non_euclidean
  )
  non_euclidean_fit <- suppressWarnings(permanova_crossed(
    non_euclidean_design, p_method = "satterthwaite"
  ))
  suppression_ok <- !non_euclidean_design$spectral$euclidean &&
    is.na(non_euclidean_fit$tests$treatment$p_value)
  rows[[length(rows) + 1L]] <- .validation_row(
    "non-Euclidean Satterthwaite p-value suppressed", suppression_ok,
    non_euclidean_design$spectral$negative_count, 1,
    non_euclidean_fit$tests$treatment$approximation_status
  )
  non_euclidean_permutation <- suppressWarnings(permanova_crossed(
    non_euclidean_design, p_method = "permutation", permutations = "all"
  ))
  non_euclidean_permutation_ok <-
    !non_euclidean_design$spectral$euclidean &&
    is.finite(non_euclidean_permutation$tests$treatment$p_value) &&
    non_euclidean_permutation$tests$treatment$p_value >= 0 &&
    non_euclidean_permutation$tests$treatment$p_value <= 1
  rows[[length(rows) + 1L]] <- .validation_row(
    "non-Euclidean crossed permutation p-value retained",
    non_euclidean_permutation_ok,
    non_euclidean_permutation$tests$treatment$p_value, 1,
    non_euclidean_permutation$tests$treatment$approximation_status
  )

  nested_data <- expand.grid(
    replicate = factor(1:3), site = factor(1:2), region = factor(1:6)
  )
  nested_data$treatment <- factor(
    ifelse(as.integer(nested_data$region) <= 3L, "control", "impact")
  )
  nested_region_effect <- c(-1.2, 0.1, 1, -0.8, 0.3, 1.3)
  nested_site_effect <- c(
    -0.4, 0.4, 0.2, -0.2, -0.1, 0.1,
    0.3, -0.3, -0.2, 0.2, 0.5, -0.5
  )
  nested_site <- interaction(
    nested_data$region, nested_data$site, drop = TRUE
  )
  nested_response <- 1.5 * (nested_data$treatment == "impact") +
    nested_region_effect[nested_data$region] +
    nested_site_effect[nested_site] + rep(c(-0.1, 0, 0.1), 12)
  nested_design <- compile_permanova_design(
    ~ treatment + (1 | region/site), nested_data, dist(nested_response)
  )
  nested_fit <- suppressWarnings(permanova_nested(
    nested_design, permutations = "all"
  ))
  nested_means <- aggregate(
    nested_response,
    list(region = nested_data$region, treatment = nested_data$treatment),
    mean
  )
  names(nested_means)[[3L]] <- "response"
  nested_reference <- unname(
    summary(aov(response ~ treatment, nested_means))[[1]][["F value"]][1]
  )
  nested_difference <- abs(
    nested_fit$tests$treatment$statistic - nested_reference
  )
  rows[[length(rows) + 1L]] <- .validation_row(
    "nested Euclidean F equals top-unit ANOVA", nested_difference < 1e-10,
    nested_difference, 1e-10, "region/site exact EMS denominator"
  )
  nested_exact_ok <- nested_fit$permutation$attainable == factorial(6) &&
    nested_fit$permutation$used == factorial(6) &&
    abs(nested_fit$tests$treatment$p_value * factorial(6) -
          round(nested_fit$tests$treatment$p_value * factorial(6))) < 1e-12
  rows[[length(rows) + 1L]] <- .validation_row(
    "exhaustive nested top-unit enumeration", nested_exact_ok,
    nested_fit$permutation$used, factorial(6), "all 6! intact hierarchies"
  )

  nested_data$site_treatment <- nested_data$site
  nested_data$observation_treatment <- nested_data$replicate
  level_response <- nested_response + rep(
    c(-0.4, 0.1, 0.3, 0.2, -0.3, 0.1), 6
  )
  level_design <- compile_permanova_design(
    ~ treatment + site_treatment + observation_treatment +
      (1 | region/site),
    nested_data, dist(level_response)
  )
  level_fit <- permanova_nested(level_design, permutations = 0)
  level_site_means <- aggregate(
    level_response,
    list(
      region = nested_data$region, site = nested_data$site,
      site_treatment = nested_data$site_treatment
    ),
    mean
  )
  names(level_site_means)[[4L]] <- "response"
  level_site_reference <- summary(aov(
    response ~ region + site_treatment, level_site_means
  ))[[1]]
  nested_site_factor <- interaction(
    nested_data$region, nested_data$site, drop = TRUE
  )
  level_observation_reference <- summary(aov(
    level_response ~ nested_site_factor + observation_treatment,
    nested_data
  ))[[1]]
  site_difference <- abs(
    level_fit$tests$site_treatment$statistic -
      level_site_reference["site_treatment", "F value"]
  )
  observation_difference <- abs(
    level_fit$tests$observation_treatment$statistic -
      level_observation_reference["observation_treatment", "F value"]
  )
  rows[[length(rows) + 1L]] <- .validation_row(
    "nested site-level F equals within-region ANOVA", site_difference < 1e-10,
    site_difference, 1e-10, "site units permute within regions"
  )
  rows[[length(rows) + 1L]] <- .validation_row(
    "nested observation-level F equals within-site ANOVA",
    observation_difference < 1e-10, observation_difference, 1e-10,
    "observations permute within nested sites"
  )

  site_exact_design <- compile_permanova_design(
    ~ site_treatment + (1 | region/site),
    nested_data, dist(level_response)
  )
  site_exact <- suppressWarnings(permanova_nested(
    site_exact_design, permutations = "all"
  ))
  site_exact_ok <- site_exact$permutation$attainable == 2 ^ 6 &&
    site_exact$permutation$used == 2 ^ 6
  rows[[length(rows) + 1L]] <- .validation_row(
    "exhaustive nested within-parent enumeration", site_exact_ok,
    site_exact$permutation$used, 2 ^ 6, "all 2^6 site arrangements"
  )

  observation_data <- expand.grid(
    replicate = factor(1:2), site = factor(1:2), region = factor(1:3)
  )
  observation_data$observation_treatment <- observation_data$replicate
  observation_response <- sin(seq_len(nrow(observation_data)))
  observation_design <- compile_permanova_design(
    ~ observation_treatment + (1 | region/site),
    observation_data, dist(observation_response)
  )
  observation_exact <- suppressWarnings(permanova_nested(
    observation_design, permutations = "all"
  ))
  observation_exact_ok <- observation_exact$permutation$attainable == 2 ^ 6 &&
    observation_exact$permutation$used == 2 ^ 6
  rows[[length(rows) + 1L]] <- .validation_row(
    "exhaustive nested within-unit enumeration", observation_exact_ok,
    observation_exact$permutation$used, 2 ^ 6,
    "all 2^6 within-site observation arrangements"
  )

  dispersion_data <- data.frame(
    treatment = factor(rep(c("control", "impact"), each = 3L))
  )
  dispersion_response <- cbind(
    c(-1.2, -0.1, 0.8, 2, 2.4, 4.1),
    c(0.1, 0.5, -0.4, 0.2, 1.3, -0.8)
  )
  dispersion_design <- compile_permanova_design(
    ~ treatment, dispersion_data, dist(dispersion_response)
  )
  dispersion_fit <- permanova_dispersion(
    dispersion_design, "treatment", type = "centroid", permutations = 0
  )
  dispersion_reference <- summary(aov(
    dispersion_fit$distances ~ dispersion_data$treatment
  ))[[1]][1L, "F value"]
  dispersion_difference <- abs(
    dispersion_fit$test$statistic - dispersion_reference
  )
  rows[[length(rows) + 1L]] <- .validation_row(
    "dispersion F equals ANOVA on distances to centroids",
    dispersion_difference < 1e-10, dispersion_difference, 1e-10,
    "independent scalar response reference"
  )

  crossed_dispersion_response <- cbind(
    crossed_response,
    sin(seq_along(crossed_response)) +
      0.4 * (crossed_data$treatment == "impact")
  )
  crossed_dispersion_design <- compile_permanova_design(
    ~ treatment + (1 | block) + (1 | treatment:block),
    crossed_data, dist(crossed_dispersion_response)
  )
  crossed_dispersion <- suppressWarnings(permanova_dispersion(
    crossed_dispersion_design, "treatment", type = "centroid",
    permutations = "all", keep_permutations = TRUE
  ))
  crossed_dispersion_structure <- .validate_crossed_design(
    .dispersion_response_design(
      crossed_dispersion_design, crossed_dispersion$distances, 1e-10
    )
  )
  crossed_dispersion_units <- .validate_crossed_permutation_design(
    .dispersion_response_design(
      crossed_dispersion_design, crossed_dispersion$distances, 1e-10
    ),
    crossed_dispersion_structure
  )
  crossed_dispersion_plan <- .crossed_block_permutation_plan(
    crossed_dispersion_units, "all", 100L
  )
  crossed_dispersion_residuals <- stats::residuals(stats::lm(
    crossed_dispersion$distances ~ crossed_data$treatment
  ))
  crossed_dispersion_direct <- vapply(
    seq_len(nrow(crossed_dispersion_plan$orders)), function(row) {
      index <- .block_permutation_index(
        crossed_dispersion_units$cell_blocks,
        crossed_dispersion_plan$orders[row, ]
      )
      reference <- summary(aov(
        crossed_dispersion_residuals[index] ~
          crossed_data$treatment * crossed_data$block
      ))[[1]]
      reference_rows <- trimws(rownames(reference))
      reference[match("crossed_data$treatment", reference_rows), "Mean Sq"] /
        reference[match(
          "crossed_data$treatment:crossed_data$block", reference_rows
        ), "Mean Sq"]
    }, numeric(1)
  )
  crossed_dispersion_difference <- max(abs(
    crossed_dispersion$test$permutation_statistics -
      crossed_dispersion_direct
  ))
  rows[[length(rows) + 1L]] <- .validation_row(
    "crossed dispersion permutations equal direct blocked residual ANOVA",
    crossed_dispersion_difference < 1e-10,
    crossed_dispersion_difference, 1e-10,
    "all 2^4 intact-cell residual arrangements"
  )

  dispersion_spectrum_ok <-
    identical(length(crossed_dispersion$eigenvalues),
              length(crossed_dispersion_design$spectral$eigenvalues)) &&
    max(abs(
      crossed_dispersion$eigenvalues -
        crossed_dispersion_design$spectral$eigenvalues
    )) < 1e-10
  rows[[length(rows) + 1L]] <- .validation_row(
    "dispersion geometry retains the complete distance spectrum",
    dispersion_spectrum_ok,
    max(abs(
      crossed_dispersion$eigenvalues -
        crossed_dispersion_design$spectral$eigenvalues
    )), 1e-10, "positive and negative axes preserved"
  )

  hybrid_data <- expand.grid(
    replicate = factor(1:2), treatment = factor(c("control", "impact")),
    site = factor(1:2), region = factor(1:3)
  )
  hybrid_region_site <- interaction(
    hybrid_data$region, hybrid_data$site, drop = TRUE
  )
  hybrid_treatment_region <- interaction(
    hybrid_data$treatment, hybrid_data$region, drop = TRUE
  )
  hybrid_treatment_site <- interaction(
    hybrid_data$treatment, hybrid_data$region, hybrid_data$site,
    drop = TRUE
  )
  hybrid_response <-
    c(-0.8, 0.1, 0.7)[hybrid_data$region] +
    seq(-0.5, 0.5, length.out = 6)[hybrid_region_site] +
    c(-0.9, 0.4, 0.7, -0.2, 0.5, -0.5)[hybrid_treatment_region] +
    seq(-0.3, 0.3, length.out = 12)[hybrid_treatment_site] +
    rep(c(-0.1, 0.1), 12)
  hybrid_design <- compile_permanova_design(
    ~ treatment + (1 | region/site) +
      (1 | treatment:region) + (1 | treatment:region:site),
    hybrid_data, dist(hybrid_response)
  )
  hybrid_fit <- suppressWarnings(permanova_hybrid(
    hybrid_design, permutations = "all", keep_permutations = TRUE
  ))
  hybrid_profile_f <- function(values) {
    means <- aggregate(
      values,
      list(treatment = hybrid_data$treatment, region = hybrid_data$region),
      mean
    )
    table <- summary(aov(x ~ treatment * region, means))[[1]]
    table_rows <- trimws(rownames(table))
    unname(table[match("treatment", table_rows), "Mean Sq"] /
             table[match("treatment:region", table_rows), "Mean Sq"])
  }
  hybrid_reference <- hybrid_profile_f(hybrid_response)
  hybrid_difference <- abs(hybrid_fit$test$statistic - hybrid_reference)
  rows[[length(rows) + 1L]] <- .validation_row(
    "hybrid treatment F equals treatment-region profile ANOVA",
    hybrid_difference < 1e-10, hybrid_difference, 1e-10,
    "exact MS[treatment:region] denominator"
  )

  hybrid_structure <- .validate_hybrid_design(hybrid_design)
  hybrid_plan <- .hybrid_permutation_plan(hybrid_structure, "all", 100L)
  hybrid_direct <- vapply(seq_len(nrow(hybrid_plan$orders)), function(row) {
    index <- .block_permutation_index(
      hybrid_structure$profile_blocks, hybrid_plan$orders[row, ]
    )
    hybrid_profile_f(hybrid_response[index])
  }, numeric(1))
  hybrid_permutation_difference <- max(abs(
    hybrid_fit$test$permutation_statistics - hybrid_direct
  ))
  rows[[length(rows) + 1L]] <- .validation_row(
    "hybrid profile permutations equal direct ANOVA",
    hybrid_permutation_difference < 1e-10,
    hybrid_permutation_difference, 1e-10,
    "all 2^3 within-region profile arrangements"
  )

  hybrid_components <- hybrid_fit$variance_components
  hybrid_component_difference <- max(abs(
    as.vector(
      hybrid_components$coefficient_matrix %*%
        hybrid_components$unrestricted
    ) - unname(hybrid_components$observed_mean_squares)
  ))
  rows[[length(rows) + 1L]] <- .validation_row(
    "hybrid variance components reconstruct selected strata",
    hybrid_component_difference < 1e-10,
    hybrid_component_difference, 1e-10,
    "five-component crossed-nested EMS system"
  )

  unbalanced_cells <- expand.grid(
    A = factor(c("a1", "a2")), B = factor(c("b1", "b2", "b3"))
  )
  unbalanced_counts <- c(2L, 4L, 3L, 2L, 5L, 3L)
  unbalanced_fixed_data <- unbalanced_cells[
    rep(seq_len(nrow(unbalanced_cells)), unbalanced_counts), , drop = FALSE
  ]
  rownames(unbalanced_fixed_data) <- NULL
  unbalanced_fixed_response <-
    1.1 * (unbalanced_fixed_data$A == "a2") +
    c(b1 = -0.5, b2 = 0.2, b3 = 0.8)[unbalanced_fixed_data$B] +
    0.2 * sin(seq_len(nrow(unbalanced_fixed_data)))
  unbalanced_fixed_design <- compile_permanova_design(
    stats::reformulate(c("A", "B")), unbalanced_fixed_data,
    dist(unbalanced_fixed_response)
  )
  unbalanced_type1 <- permanova_unbalanced(
    unbalanced_fixed_design, hypothesis = "I", permutations = 0
  )
  unbalanced_reference_data <- unbalanced_fixed_data
  unbalanced_reference_data$response <- unbalanced_fixed_response
  unbalanced_type1_reference <- summary(aov(
    stats::reformulate(c("A", "B"), response = "response"),
    unbalanced_reference_data
  ))[[1]]
  unbalanced_type1_difference <- max(abs(
    vapply(unbalanced_type1$tests, `[[`, numeric(1), "statistic") -
      unbalanced_type1_reference[c("A", "B"), "F value"]
  ))
  rows[[length(rows) + 1L]] <- .validation_row(
    "unbalanced Type I statistics equal sequential ANOVA",
    unbalanced_type1_difference < 1e-10,
    unbalanced_type1_difference, 1e-10,
    "unequal complete two-factor cells"
  )

  unbalanced_type2_forward <- permanova_unbalanced(
    unbalanced_fixed_design, hypothesis = "II", permutations = 0
  )
  unbalanced_type2_reverse <- permanova_unbalanced(
    compile_permanova_design(
      stats::reformulate(c("B", "A")), unbalanced_fixed_data,
      dist(unbalanced_fixed_response)
    ),
    hypothesis = "II", permutations = 0
  )
  unbalanced_type2_difference <- max(abs(c(
    unbalanced_type2_forward$tests$A$statistic -
      unbalanced_type2_reverse$tests$A$statistic,
    unbalanced_type2_forward$tests$B$statistic -
      unbalanced_type2_reverse$tests$B$statistic
  )))
  rows[[length(rows) + 1L]] <- .validation_row(
    "unbalanced additive Type II statistics are order invariant",
    unbalanced_type2_difference < 1e-10,
    unbalanced_type2_difference, 1e-10,
    "full-versus-drop marginal hypotheses"
  )

  unbalanced_type3_data <- unbalanced_fixed_data[
    unbalanced_fixed_data$B != "b3", , drop = FALSE
  ]
  unbalanced_type3_data <- droplevels(unbalanced_type3_data)
  attr(unbalanced_type3_data$A, "contrasts") <- stats::contr.sum(2)
  attr(unbalanced_type3_data$B, "contrasts") <- stats::contr.sum(2)
  unbalanced_type3_response <- unbalanced_fixed_response[
    unbalanced_fixed_data$B != "b3"
  ] + 0.6 * (
    unbalanced_type3_data$A == "a2" &
      unbalanced_type3_data$B == "b2"
  )
  unbalanced_type3_design <- compile_permanova_design(
    stats::reformulate(c("A", "B", "A:B")),
    unbalanced_type3_data, dist(unbalanced_type3_response)
  )
  unbalanced_type3 <- permanova_unbalanced(
    unbalanced_type3_design, hypothesis = "III", permutations = 0
  )
  type3_X <- unbalanced_type3_design$X
  type3_assignment <- attr(type3_X, "assign")
  type3_full <- stats::lm.fit(type3_X, unbalanced_type3_response)
  type3_full_rss <- sum(type3_full$residuals ^ 2)
  type3_full_df <- nrow(type3_X) - type3_full$rank
  type3_direct <- vapply(
    seq_along(unbalanced_type3_design$fixed_effects), function(index) {
      reduced <- stats::lm.fit(
        type3_X[, type3_assignment != index, drop = FALSE],
        unbalanced_type3_response
      )
      numerator_df <- type3_full$rank - reduced$rank
      ((sum(reduced$residuals ^ 2) - type3_full_rss) / numerator_df) /
        (type3_full_rss / type3_full_df)
    }, numeric(1)
  )
  unbalanced_type3_difference <- max(abs(
    vapply(unbalanced_type3$tests, `[[`, numeric(1), "statistic") -
      type3_direct
  ))
  rows[[length(rows) + 1L]] <- .validation_row(
    "unbalanced Type III statistics equal full-versus-drop references",
    unbalanced_type3_difference < 1e-10,
    unbalanced_type3_difference, 1e-10,
    "explicit sum-to-zero contrasts"
  )

  unbalanced_exact_data <- data.frame(
    group = factor(c("a", "a", "b", "b", "b")),
    response = c(-1.1, 0.3, 0.4, 1.8, 3.2)
  )
  unbalanced_exact_design <- compile_permanova_design(
    ~ group, unbalanced_exact_data, dist(unbalanced_exact_data$response)
  )
  unbalanced_exact <- suppressWarnings(permanova_unbalanced(
    unbalanced_exact_design, hypothesis = "III", permutations = "all",
    keep_permutations = TRUE, max_exact = 1000L
  ))
  unbalanced_orders <- .all_orders(nrow(unbalanced_exact_data))
  unbalanced_reduced <- stats::lm(response ~ 1, unbalanced_exact_data)
  unbalanced_direct <- apply(unbalanced_orders, 1L, function(index) {
    permuted_response <- stats::fitted(unbalanced_reduced) +
      stats::residuals(unbalanced_reduced)[index]
    summary(aov(
      permuted_response ~ group, unbalanced_exact_data
    ))[[1]][1L, "F value"]
  })
  unbalanced_permutation_difference <- max(abs(
    unbalanced_exact$tests$group$permutation_statistics -
      unbalanced_direct
  ))
  rows[[length(rows) + 1L]] <- .validation_row(
    "unbalanced exact residual permutations equal direct ANOVA",
    unbalanced_permutation_difference < 1e-10,
    unbalanced_permutation_difference, 1e-10,
    "all 5! Freedman-Lane arrangements"
  )

  incomplete_rows <- !(
    unbalanced_fixed_data$A == "a2" &
      unbalanced_fixed_data$B == "b3"
  )
  incomplete_data <- unbalanced_fixed_data[
    incomplete_rows, , drop = FALSE
  ]
  incomplete_data <- droplevels(incomplete_data)
  attr(incomplete_data$A, "contrasts") <- stats::contr.sum(2)
  attr(incomplete_data$B, "contrasts") <- stats::contr.sum(3)
  incomplete_design <- compile_permanova_design(
    stats::reformulate(c("A", "B", "A:B")), incomplete_data,
    dist(unbalanced_fixed_response[
      !(unbalanced_fixed_data$A == "a2" & unbalanced_fixed_data$B == "b3")
    ])
  )
  incomplete_message <- tryCatch({
    permanova_unbalanced(
      incomplete_design, hypothesis = "III", permutations = 0
    )
    "not rejected"
  }, error = conditionMessage)
  rows[[length(rows) + 1L]] <- .validation_row(
    "aliased incomplete Type III hypothesis rejects explicitly",
    grepl("full-rank fixed model", incomplete_message),
    as.numeric(grepl("full-rank fixed model", incomplete_message)), 1,
    incomplete_message
  )

  unequal_sizes <- c(2L, 2L, 3L, 3L, 4L, 4L)
  unequal_mixed_data <- data.frame(
    site = factor(rep(seq_along(unequal_sizes), unequal_sizes)),
    treatment = factor(rep(rep(c("control", "impact"), 3L), unequal_sizes))
  )
  unequal_unit_effect <- c(-0.8, 0.4, -0.2, 0.7, -0.5, 1.0)
  unequal_mixed_response <- 0.6 * (unequal_mixed_data$treatment == "impact") +
    unequal_unit_effect[unequal_mixed_data$site] +
    stats::ave(seq_len(nrow(unequal_mixed_data)), unequal_mixed_data$site,
               FUN = function(x) seq(-0.15, 0.15, length.out = length(x)))
  unequal_mixed_design <- compile_permanova_design(
    ~ treatment + (1 | site), unequal_mixed_data, dist(unequal_mixed_response)
  )
  unequal_mixed_fit <- suppressWarnings(permanova_unequal_random_intercept(
    unequal_mixed_design, permutations = "all", keep_permutations = TRUE
  ))
  unequal_means <- aggregate(
    unequal_mixed_response,
    list(site = unequal_mixed_data$site,
         treatment = unequal_mixed_data$treatment), mean
  )
  unequal_means <- unequal_means[
    order(as.integer(unequal_means$site)), , drop = FALSE
  ]
  unequal_reference <- summary(aov(x ~ treatment, unequal_means))[[1]][1L, "F value"]
  unequal_difference <- abs(unequal_mixed_fit$test$statistic - unequal_reference)
  rows[[length(rows) + 1L]] <- .validation_row(
    "unequal random-intercept F equals cluster-centroid ANOVA",
    unequal_difference < 1e-10, unequal_difference, 1e-10,
    "cluster-weighted estimand"
  )
  unequal_structure <- .validate_unequal_random_intercept(unequal_mixed_design)
  unequal_plan <- .stratified_permutation_plan(
    unequal_structure$size_strata, "all", 100L
  )
  unequal_reduced <- stats::lm(x ~ 1, unequal_means)
  unequal_direct <- apply(unequal_plan$orders, 1L, function(index) {
    response <- stats::fitted(unequal_reduced) +
      stats::residuals(unequal_reduced)[index]
    summary(aov(response ~ treatment, unequal_means))[[1]][1L, "F value"]
  })
  unequal_permutation_difference <- max(abs(
    unequal_mixed_fit$test$permutation_statistics - unequal_direct
  ))
  rows[[length(rows) + 1L]] <- .validation_row(
    "unequal random-intercept permutations equal direct centroid ANOVA",
    unequal_permutation_difference < 1e-10,
    unequal_permutation_difference, 1e-10,
    "all 2! x 2! x 2! size-stratified arrangements"
  )

  calibration_smoke <- audit_unequal_random_intercept_calibration(
    batches = 1L, simulations_per_batch = 2L, permutations = 9L,
    engine_checks = 2L, seed = seed + 20260000L
  )
  rows[[length(rows) + 1L]] <- .validation_row(
    "unique assignment audit equals exhaustive distance engine",
    calibration_smoke$max_engine_exact_p_difference < 1e-12,
    calibration_smoke$max_engine_exact_p_difference, 1e-12,
    "216 unique labels equal 24^3 profile permutations"
  )
  phase14_strata <- factor(rep(seq_len(3L), each = 4L))
  phase14_plan <- .stratified_permutation_plan(
    phase14_strata, 5L, 1000L
  )
  phase14_H0 <- matrix(1 / 12, 12, 12)
  phase14_R0 <- diag(12) - phase14_H0
  phase14_identity_error <- max(apply(
    phase14_plan$orders, 1L, function(index) {
      P <- diag(12)[index, , drop = FALSE]
      max(abs(
        phase14_H0 + phase14_R0[index, , drop = FALSE] - P
      ))
    }
  ))
  rows[[length(rows) + 1L]] <- .validation_row(
    "intercept-only Freedman-Lane equals raw centroid permutation",
    phase14_identity_error < 1e-12,
    phase14_identity_error, 1e-12,
    "algebraic exchangeability identity"
  )
  phase14_seed_spacing <- if (
    length(calibration_smoke$seed_design$batch_seeds) > 1L
  ) min(diff(calibration_smoke$seed_design$batch_seeds)) else 1000000L
  rows[[length(rows) + 1L]] <- .validation_row(
    "calibration batches use disjoint deterministic seed ranges",
    phase14_seed_spacing == 1000000L,
    phase14_seed_spacing, 1000000L,
    calibration_smoke$seed_design$note
  )

  phase12_nested_source <- .phase12_nested_data()
  phase12_nested_data <- droplevels(
    phase12_nested_source[as.integer(phase12_nested_source$region) <= 6L, ]
  )
  phase12_nested_response <-
    c(-0.9, 0.2, -0.4, 0.8, -0.2, 1.1)[phase12_nested_data$region] +
    seq(-0.35, 0.35, length.out = nlevels(phase12_nested_data$site))[
      phase12_nested_data$site
    ] + sin(seq_len(nrow(phase12_nested_data))) / 20
  phase12_nested_design <- compile_permanova_design(
    ~ treatment + (1 | region/site), phase12_nested_data,
    dist(phase12_nested_response)
  )
  phase12_nested_fit <- suppressWarnings(permanova_unbalanced_nested(
    phase12_nested_design, permutations = "all", keep_permutations = TRUE
  ))
  phase12_top_means <- aggregate(
    phase12_nested_response,
    list(region = phase12_nested_data$region,
         treatment = phase12_nested_data$treatment), mean
  )
  phase12_top_means <- phase12_top_means[
    order(as.integer(phase12_top_means$region)), , drop = FALSE
  ]
  names(phase12_top_means)[[3L]] <- "response"
  phase12_nested_reference <- summary(aov(
    response ~ treatment, phase12_top_means
  ))[[1]][1L, "F value"]
  phase12_nested_difference <- abs(
    phase12_nested_fit$test$statistic - phase12_nested_reference
  )
  rows[[length(rows) + 1L]] <- .validation_row(
    "unbalanced nested F equals top-centroid ANOVA",
    phase12_nested_difference < 1e-10,
    phase12_nested_difference, 1e-10, "topology-conditional top units"
  )
  phase12_nested_structure <- .validate_unbalanced_nested(
    phase12_nested_design
  )
  phase12_nested_plan <- .stratified_permutation_plan(
    phase12_nested_structure$topology_strata, "all", 100L
  )
  phase12_nested_reduced <- stats::lm(response ~ 1, phase12_top_means)
  phase12_nested_direct <- apply(
    phase12_nested_plan$orders, 1L, function(index) {
      permutation_data <- phase12_top_means
      permutation_data$response <- stats::fitted(phase12_nested_reduced) +
        stats::residuals(phase12_nested_reduced)[index]
      summary(aov(response ~ treatment, permutation_data))[[1]][1L, "F value"]
    }
  )
  phase12_nested_permutation_difference <- max(abs(
    phase12_nested_fit$test$permutation_statistics - phase12_nested_direct
  ))
  rows[[length(rows) + 1L]] <- .validation_row(
    "unbalanced nested permutations equal direct top-centroid ANOVA",
    phase12_nested_permutation_difference < 1e-10,
    phase12_nested_permutation_difference, 1e-10,
    "all paired-topology arrangements"
  )

  phase12_block_source <- .phase12_block_data()
  phase12_block_data <- droplevels(
    phase12_block_source[as.integer(phase12_block_source$block) <= 4L, ]
  )
  phase12_block_cell <- interaction(
    phase12_block_data$treatment, phase12_block_data$block, drop = TRUE
  )
  phase12_block_response <-
    c(-0.6, 0.1, 0.4, -0.2)[phase12_block_data$block] +
    c(A = 0, B = 0.4, C = 0.8)[phase12_block_data$treatment] +
    sin(as.integer(phase12_block_cell) * 1.7)[phase12_block_cell] +
    cos(seq_len(nrow(phase12_block_data))) / 20
  phase12_block_design <- compile_permanova_design(
    ~ treatment + (1 | block) + (1 | treatment:block),
    phase12_block_data, dist(phase12_block_response)
  )
  phase12_block_fit <- suppressWarnings(permanova_unbalanced_block(
    phase12_block_design, permutations = "all", keep_permutations = TRUE
  ))
  phase12_block_structure <- .validate_unbalanced_block(phase12_block_design)
  phase12_cell_means <- aggregate(
    phase12_block_response,
    list(treatment = phase12_block_data$treatment,
         block = phase12_block_data$block), mean
  )
  phase12_cell_means <- phase12_cell_means[match(
    paste(phase12_block_structure$cell_mapping$treatment,
          phase12_block_structure$cell_mapping$block),
    paste(phase12_cell_means$treatment, phase12_cell_means$block)
  ), , drop = FALSE]
  names(phase12_cell_means)[[3L]] <- "response"
  phase12_block_table <- summary(aov(
    response ~ block + treatment, phase12_cell_means
  ))[[1]]
  phase12_block_reference <- phase12_block_table["treatment", "Mean Sq"] /
    phase12_block_table["Residuals", "Mean Sq"]
  phase12_block_difference <- abs(
    phase12_block_fit$test$statistic - phase12_block_reference
  )
  rows[[length(rows) + 1L]] <- .validation_row(
    "unbalanced block F equals cell-centroid blocked ANOVA",
    phase12_block_difference < 1e-10,
    phase12_block_difference, 1e-10, "incomplete cell incidence"
  )
  phase12_block_plan <- .stratified_permutation_plan(
    phase12_block_structure$exchangeability_strata, "all", 1000L
  )
  phase12_block_direct <- apply(
    phase12_block_plan$orders, 1L, function(index) {
      permutation_data <- phase12_cell_means
      permutation_data$response <- phase12_cell_means$response[index]
      table <- summary(aov(
        response ~ block + treatment, permutation_data
      ))[[1]]
      table["treatment", "Mean Sq"] / table["Residuals", "Mean Sq"]
    }
  )
  phase12_block_permutation_difference <- max(abs(
    phase12_block_fit$test$permutation_statistics - phase12_block_direct
  ))
  rows[[length(rows) + 1L]] <- .validation_row(
    "unbalanced block permutations equal direct cell-centroid ANOVA",
    phase12_block_permutation_difference < 1e-10,
    phase12_block_permutation_difference, 1e-10,
    "all within-block intact-cell arrangements"
  )

  contrast_data <- data.frame(
    group = factor(c("A", "A", "B", "B", "C"),
                   levels = c("A", "B", "C"))
  )
  contrast_response <- c(-1.2, 0.2, 0.4, 1.8, 2.5)
  contrast_matrix <- rbind(
    `A - B` = c(A = 1, B = -1, C = 0),
    `AB - C` = c(A = 0.5, B = 0.5, C = -1)
  )
  contrast_design <- compile_permanova_design(
    ~ group, contrast_data, dist(contrast_response)
  )
  contrast_fit <- suppressWarnings(permanova_contrasts(
    contrast_design, contrast_matrix, permutations = "all",
    keep_permutations = TRUE
  ))
  contrast_means <- tapply(
    contrast_response, contrast_data$group, mean
  )
  contrast_counts <- table(contrast_data$group)
  contrast_full <- stats::lm(contrast_response ~ group, contrast_data)
  contrast_residual_ms <- sum(stats::residuals(contrast_full) ^ 2) /
    stats::df.residual(contrast_full)
  contrast_direct_f <- apply(contrast_matrix, 1L, function(coefficients) {
    (sum(coefficients * contrast_means) ^ 2 /
       sum(coefficients ^ 2 / contrast_counts)) / contrast_residual_ms
  })
  contrast_f_difference <- max(abs(
    contrast_fit$table$pseudo_f - contrast_direct_f
  ))
  rows[[length(rows) + 1L]] <- .validation_row(
    "planned contrast F equals scalar contrast formula",
    contrast_f_difference < 1e-10,
    contrast_f_difference, 1e-10, "unequal one-factor group counts"
  )
  contrast_coefficients <- contrast_matrix[1L, ]
  contrast_X <- stats::model.matrix(~ 0 + group, contrast_data)
  contrast_Q <- qr.Q(
    qr(matrix(contrast_coefficients, ncol = 1L)), complete = TRUE
  )
  contrast_reduced_X <- contrast_X %*% contrast_Q[, -1L, drop = FALSE]
  contrast_reduced <- stats::lm.fit(
    contrast_reduced_X, contrast_response
  )
  contrast_fitted <- as.vector(
    contrast_reduced_X %*% contrast_reduced$coefficients
  )
  contrast_residual <- contrast_response - contrast_fitted
  contrast_orders <- .all_orders(nrow(contrast_data))
  contrast_direct_permutations <- apply(
    contrast_orders, 1L, function(index) {
      response <- contrast_fitted + contrast_residual[index]
      means <- tapply(response, contrast_data$group, mean)
      ss <- sum(contrast_coefficients * means) ^ 2 /
        sum(contrast_coefficients ^ 2 / contrast_counts)
      full <- stats::lm(response ~ group, contrast_data)
      ss / (sum(stats::residuals(full) ^ 2) / stats::df.residual(full))
    }
  )
  contrast_permutation_difference <- max(abs(
    contrast_fit$tests[[1L]]$permutation_statistics -
      contrast_direct_permutations
  ))
  rows[[length(rows) + 1L]] <- .validation_row(
    "planned contrast permutations equal constrained scalar fit",
    contrast_permutation_difference < 1e-10,
    contrast_permutation_difference, 1e-10,
    "all 5! contrast-specific Freedman-Lane arrangements"
  )
  contrast_adjustment_difference <- max(abs(
    contrast_fit$table$p_adjusted -
      stats::p.adjust(contrast_fit$table$p_perm, method = "holm")
  ))
  rows[[length(rows) + 1L]] <- .validation_row(
    "contrast family adjustment equals stats::p.adjust",
    contrast_adjustment_difference < 1e-12,
    contrast_adjustment_difference, 1e-12,
    "Holm adjustment under arbitrary contrast dependence"
  )

  marginal_validation <- simulate_marginal_contrast_validation(
    simulations = 2L, permutations = 7L, seed = seed + 20260719L
  )
  rows[[length(rows) + 1L]] <- .validation_row(
    "multifactor marginal contrast F equals scalar GLH",
    marginal_validation$max_engine_direct_f_difference < 1e-10,
    marginal_validation$max_engine_direct_f_difference, 1e-10,
    "equal nuisance-grid weighting"
  )
  rows[[length(rows) + 1L]] <- .validation_row(
    "multifactor marginal permutations equal response reconstruction",
    marginal_validation$max_engine_direct_permutation_statistic_difference < 1e-10,
    marginal_validation$max_engine_direct_permutation_statistic_difference,
    1e-10, "contrast-specific Freedman-Lane"
  )

  mixed_contrast_data <- expand.grid(
    replicate = factor(1:2), treatment = factor(c("A", "B", "C")),
    block = factor(1:4)
  )
  mixed_contrast_cell <- interaction(
    mixed_contrast_data$treatment, mixed_contrast_data$block, drop = TRUE
  )
  mixed_contrast_response <-
    c(-1, 0, 1, 0.2, 1, 2, 1, 2, 4, -0.5, 1, 3)[
      mixed_contrast_cell
    ] + c(-0.1, 0.1)
  mixed_contrast_design <- compile_permanova_design(
    ~ treatment + (1 | block) + (1 | treatment:block),
    mixed_contrast_data, dist(mixed_contrast_response)
  )
  mixed_contrast_fit <- permanova_mixed_contrasts(
    mixed_contrast_design, c(A = 1, B = -1, C = 0),
    permutations = 0
  )
  mixed_contrast_means <- aggregate(
    mixed_contrast_response,
    list(treatment = mixed_contrast_data$treatment,
         block = mixed_contrast_data$block), mean
  )
  mixed_contrast_reference <- .scalar_mixed_contrast_reference(
    mixed_contrast_means$x, mixed_contrast_means$treatment,
    c(A = 1, B = -1, C = 0), mixed_contrast_means$block
  )
  mixed_contrast_difference <- abs(
    mixed_contrast_fit$table$pseudo_f - mixed_contrast_reference$F
  )
  rows[[length(rows) + 1L]] <- .validation_row(
    "mixed contrast F equals cell-centroid scalar GLH",
    mixed_contrast_difference < 1e-10,
    mixed_contrast_difference, 1e-10,
    "declared treatment-by-block denominator"
  )

  maxT_fit <- permanova_maxT(
    contrast_design, contrast_matrix, permutations = "all",
    keep_permutations = TRUE
  )
  maxT_structure <- .validate_contrast_design(contrast_design)
  maxT_plan <- .permutation_orders(contrast_design$n, "all", 1000L)
  maxT_projectors <- lapply(seq_len(nrow(contrast_matrix)), function(index) {
    .contrast_projector(
      contrast_matrix[index, ], maxT_structure$group,
      maxT_structure$counts, maxT_structure$H_full, 1e-10
    )$H
  })
  maxT_direct <- t(apply(maxT_plan$orders, 1L, function(index) {
    G <- contrast_design$G[index, index, drop = FALSE]
    denominator <- trace_hgh(maxT_structure$H_denominator, G) /
      maxT_structure$denominator_df
    vapply(maxT_projectors, function(H) {
      trace_hgh(H, G) / denominator
    }, numeric(1))
  }))
  maxT_difference <- max(abs(
    unname(maxT_fit$permutation$statistic_matrix) - maxT_direct
  ))
  rows[[length(rows) + 1L]] <- .validation_row(
    "global-null max-F joint distribution equals direct enumeration",
    maxT_difference < 1e-10, maxT_difference, 1e-10,
    "all 5! shared raw-observation arrangements"
  )

  span_treatment <- matrix(c(0, 0, 0, 1, 1, 1), ncol = 1L)
  span_strata <- factor(c(2, 2, 2, 3, 3, 3))
  span_rank <- .stratified_treatment_rank(
    span_treatment, span_strata, 1e-10
  )
  rows[[length(rows) + 1L]] <- .validation_row(
    "stratified span guard rejects numerical centering dust",
    as.integer(span_rank) == 0L,
    attr(span_rank, "projected_singular_values")[[1L]],
    attr(span_rank, "absolute_floor"),
    "absolute SVD floor relative to original treatment design"
  )

  phase19_validation <- simulate_unbalanced_mixed_contrast_validation(
    simulations = 2L, permutations = 7L, seed = seed + 20260720L
  )
  rows[[length(rows) + 1L]] <- .validation_row(
    "unbalanced mixed contrast F equals cluster-centroid GLH",
    phase19_validation$max_engine_classical_f_difference < 1e-10,
    phase19_validation$max_engine_classical_f_difference, 1e-10,
    "exact-size-stratified unit-centroid estimand"
  )
  rows[[length(rows) + 1L]] <- .validation_row(
    "unbalanced mixed contrast permutations equal reconstruction",
    phase19_validation$max_engine_direct_permutation_difference < 1e-10,
    phase19_validation$max_engine_direct_permutation_difference, 1e-10,
    "contrast-specific restricted Freedman-Lane"
  )

  phase20_validation <- simulate_adjusted_contrast_validation(
    simulations = 2L, permutations = 7L, seed = seed + 20260721L
  )
  rows[[length(rows) + 1L]] <- .validation_row(
    "covariate-adjusted contrast F equals scalar GLH",
    phase20_validation$max_engine_classical_f_difference < 1e-10,
    phase20_validation$max_engine_classical_f_difference, 1e-10,
    "explicit continuous reference with categorical nuisance weighting"
  )
  rows[[length(rows) + 1L]] <- .validation_row(
    "covariate-adjusted permutations equal response reconstruction",
    phase20_validation$max_engine_direct_permutation_difference < 1e-10,
    phase20_validation$max_engine_direct_permutation_difference, 1e-10,
    "contrast-specific unrestricted Freedman-Lane"
  )

  condition_treatment <- factor(rep(c("A", "B", "C"), each = 10L))
  condition_covariate <- as.numeric(condition_treatment) +
    1e-9 * sin(seq_along(condition_treatment))
  condition_response <- cos(seq_along(condition_treatment))
  condition_design <- compile_permanova_design(
    ~ condition_treatment + condition_covariate,
    data.frame(condition_treatment, condition_covariate, condition_response),
    stats::dist(condition_response)
  )
  condition_message <- tryCatch({
    permanova_adjusted_contrasts(
      condition_design, "condition_treatment",
      c(A = 1, B = -1, C = 0),
      at = c(condition_covariate = 2), permutations = 0
    )
    ""
  }, error = conditionMessage)
  condition_guard <- grepl("ill-conditioned", condition_message, fixed = TRUE)
  rows[[length(rows) + 1L]] <- .validation_row(
    "covariate-adjusted near-alias rejects informatively",
    condition_guard, as.numeric(condition_guard), 1,
    condition_message
  )

  phase21_validation <- simulate_mixed_maxT_validation(
    simulations = 2L, permutations = 7L, seed = seed + 20260722L
  )
  rows[[length(rows) + 1L]] <- .validation_row(
    "mixed global-null max-F observed statistics equal scalar unit GLH",
    phase21_validation$max_engine_scalar_f_difference < 1e-10,
    phase21_validation$max_engine_scalar_f_difference, 1e-10,
    "exact-size-stratified experimental-unit estimand"
  )
  rows[[length(rows) + 1L]] <- .validation_row(
    "mixed global-null max-F joint distribution equals scalar reconstruction",
    phase21_validation$max_engine_scalar_statistic_difference < 1e-10,
    phase21_validation$max_engine_scalar_statistic_difference, 1e-10,
    "shared raw experimental-unit permutation"
  )

  phase25_validation <- simulate_unbalanced_nested_lower_validation(
    simulations = 2L, permutations = 7L, seed = seed + 20260723L
  )
  for (index in seq_len(nrow(phase25_validation$agreement))) {
    agreement <- phase25_validation$agreement[index, ]
    rows[[length(rows) + 1L]] <- .validation_row(
      paste0("unbalanced nested ", agreement$structure,
             "-level F equals scalar GLH"),
      agreement$max_f_difference < 1e-10,
      agreement$max_f_difference, 1e-10,
      "lower-level nuisance-adjusted unit model"
    )
    rows[[length(rows) + 1L]] <- .validation_row(
      paste0("unbalanced nested ", agreement$structure,
             "-level permutations equal reconstruction"),
      agreement$max_permutation_difference < 1e-10,
      agreement$max_permutation_difference, 1e-10,
      "restricted lower-level Freedman-Lane"
    )
  }

  phase26_validation <- simulate_unbalanced_hybrid_validation(
    simulations = 2L, permutations = 7L, seed = seed + 20260724L
  )
  rows[[length(rows) + 1L]] <- .validation_row(
    "unbalanced hybrid F equals equal-site profile GLH",
    phase26_validation$max_engine_classical_f_difference < 1e-10,
    phase26_validation$max_engine_classical_f_difference, 1e-10,
    "equal-site-weight profile centroids with additive region adjustment"
  )
  rows[[length(rows) + 1L]] <- .validation_row(
    "unbalanced hybrid permutations equal profile reconstruction",
    phase26_validation$max_engine_direct_permutation_difference < 1e-10,
    phase26_validation$max_engine_direct_permutation_difference, 1e-10,
    "region-by-exact-topology raw profile permutation"
  )

  phase27_validation <- simulate_random_slope_validation(
    simulations = 2L, permutations = 7L, seed = seed + 20260725L
  )
  rows[[length(rows) + 1L]] <- .validation_row(
    "random-slope F equals scalar paired-difference test",
    phase27_validation$max_engine_scalar_f_difference < 1e-10,
    phase27_validation$max_engine_scalar_f_difference, 1e-10,
    "complete-subject condition differences"
  )
  rows[[length(rows) + 1L]] <- .validation_row(
    "random-slope sign flips equal direct reconstruction",
    phase27_validation$max_engine_direct_permutation_difference < 1e-10,
    phase27_validation$max_engine_direct_permutation_difference, 1e-10,
    "whole-subject central-symmetry sign flips"
  )
  rows[[length(rows) + 1L]] <- .validation_row(
    "random-slope trace covariance equals scalar method of moments",
    phase27_validation$max_engine_scalar_covariance_difference < 1e-10,
    phase27_validation$max_engine_scalar_covariance_difference, 1e-10,
    "residual-corrected intercept-slope covariance"
  )

  phase28_validation <- simulate_strong_pairwise_validation(
    simulations = 2L, permutations = 7L, seed = seed + 20260726L
  )
  rows[[length(rows) + 1L]] <- .validation_row(
    "strong closed pairwise intersection statistics equal reconstruction",
    phase28_validation$max_direct_statistic_difference < 1e-10,
    phase28_validation$max_direct_statistic_difference, 1e-10,
    "intersection-specific maximum pseudo-F"
  )
  rows[[length(rows) + 1L]] <- .validation_row(
    "strong closed pairwise p-values equal reconstruction",
    phase28_validation$max_direct_p_difference < 1e-12,
    phase28_validation$max_direct_p_difference, 1e-12,
    "closed testing over separately valid intersections"
  )
  rows[[length(rows) + 1L]] <- .validation_row(
    "partial-null pairwise permutations keep non-null groups fixed",
    phase28_validation$component_preservation_violations == 0L,
    phase28_validation$component_preservation_violations, 0,
    "connected equality-component label reassignment"
  )

  capabilities <- permanova_capabilities()
  capability_engines_resolve <- all(vapply(
    capabilities$engine, exists, logical(1), mode = "function",
    envir = asNamespace("distmer"), inherits = FALSE
  ))
  rows[[length(rows) + 1L]] <- .validation_row(
    "capability registry resolves every primary engine",
    nrow(capabilities) == 21L && !anyDuplicated(capabilities$engine) &&
      capability_engines_resolve,
    nrow(capabilities), 21,
    "unique public primary-engine registry"
  )
  capability_contracts_complete <- all(nzchar(capabilities$denominator)) &&
    all(nzchar(capabilities$permutation_contract)) &&
    all(nzchar(capabilities$explicit_exclusions)) &&
    grepl("weak FWER", capabilities$error_control[
      capabilities$engine == "permanova_maxT"
    ]) &&
    grepl("strong FWER", capabilities$error_control[
      capabilities$engine == "permanova_pairwise_strong"
    ]) &&
    grepl("weak FWER", capabilities$error_control[
      capabilities$engine == "permanova_mixed_maxT"
    ])
  rows[[length(rows) + 1L]] <- .validation_row(
    "capability registry states denominator permutation and exclusions",
    capability_contracts_complete,
    sum(!nzchar(capabilities$explicit_exclusions)), 0,
    "shared max-F remains weak; closed pairwise row states strong control"
  )

  if (include_vegan && requireNamespace("vegan", quietly = TRUE)) {
    vegan_fit <- vegan::adonis2(stats::dist(y) ~ group, permutations = 0)
    vegan_difference <- abs(one$table$pseudo_f[[1]] - vegan_fit$F[[1L]])
    rows[[length(rows) + 1L]] <- .validation_row(
      "one-way F equals vegan::adonis2", vegan_difference < 1e-10,
      vegan_difference, 1e-10, "numerical comparison"
    )
    vegan_dispersion <- vegan::betadisper(
      dist(dispersion_response), dispersion_data$treatment,
      type = "centroid"
    )
    vegan_dispersion_difference <- max(abs(
      dispersion_fit$distances - vegan_dispersion$distances
    ))
    rows[[length(rows) + 1L]] <- .validation_row(
      "dispersion distances equal vegan::betadisper",
      vegan_dispersion_difference < 1e-10,
      vegan_dispersion_difference, 1e-10,
      "centroid geometry comparison"
    )
  } else {
    rows[[length(rows) + 1L]] <- .validation_row(
      "one-way F equals vegan::adonis2", TRUE, NA_real_, NA_real_,
      "skipped: vegan unavailable or comparison disabled"
    )
    rows[[length(rows) + 1L]] <- .validation_row(
      "dispersion distances equal vegan::betadisper", TRUE,
      NA_real_, NA_real_,
      "skipped: vegan unavailable or comparison disabled"
    )
  }

  report <- do.call(rbind, rows)
  rownames(report) <- NULL
  class(report) <- c("permanova_validation", class(report))
  report
}

#' @export
print.permanova_validation <- function(x, ...) {
  NextMethod("print", x, row.names = FALSE)
  cat("Overall:", if (all(x$passed)) "PASS" else "FAIL", "\n")
  invisible(x)
}

.scalar_mixed_contrast_reference <- function(y, treatment, contrast,
                                             strata = NULL) {
  data <- data.frame(treatment = treatment)
  if (!is.null(strata)) data$strata <- strata
  X <- if (is.null(strata)) {
    stats::model.matrix(~ treatment, data)
  } else {
    stats::model.matrix(~ strata + treatment, data)
  }
  beta <- qr.solve(X, y)
  mean_rows <- do.call(rbind, lapply(levels(treatment), function(level) {
    colMeans(X[treatment == level, , drop = FALSE])
  }))
  L <- as.numeric(contrast %*% mean_rows)
  estimate <- sum(L * beta)
  numerator <- estimate ^ 2 /
    as.numeric(L %*% solve(crossprod(X), L))
  residual <- y - X %*% beta
  denominator_df <- length(y) - qr(X)$rank
  denominator <- sum(residual ^ 2) / denominator_df
  list(F = numerator / denominator, df = denominator_df)
}

#' Validate mixed-design planned contrasts
#'
#' Simulates a three-treatment randomized-block experiment with a genuine
#' treatment-by-block component. The planned A-minus-B contrast is checked
#' against a scalar general-linear-hypothesis test on treatment-by-block cell
#' means under both the null and a location alternative.
#'
#' @param simulations Number of simulations per scenario.
#' @param permutations Permutations per fit.
#' @param blocks Number of random blocks.
#' @param replicates Observations per treatment-by-block cell.
#' @param effect A-minus-B location shift in the power scenario.
#' @param alpha Rejection threshold.
#' @param seed Reproducibility seed.
#' @param keep_results Retain simulation-level results.
#' @return A mixed-contrast validation report.
#' @export
simulate_mixed_contrast_validation <- function(
    simulations = 500L, permutations = 99L, blocks = 8L,
    replicates = 2L, effect = 0.8, alpha = 0.05,
    seed = 20260719L, keep_results = FALSE) {
  simulations <- .positive_integer(simulations, "simulations")
  blocks <- .positive_integer(blocks, "blocks")
  replicates <- .positive_integer(replicates, "replicates")
  if (blocks < 3L || replicates < 2L) {
    stop("At least three blocks and two observations per cell are required.",
         call. = FALSE)
  }
  if (length(permutations) != 1L || !is.finite(permutations) ||
      permutations < 0 || permutations != as.integer(permutations)) {
    stop("`permutations` must be one non-negative integer.", call. = FALSE)
  }
  if (!is.finite(effect) || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`effect` must be finite and `alpha` must lie between zero and one.",
         call. = FALSE)
  }
  template <- expand.grid(
    replicate = factor(seq_len(replicates)),
    treatment = factor(c("A", "B", "C"), levels = c("A", "B", "C")),
    block = factor(seq_len(blocks))
  )
  contrast <- c(A = 1, B = -1, C = 0)
  scenarios <- c(null = 0, location_shift = effect)
  rows <- vector("list", length(scenarios) * simulations)
  position <- 0L
  start <- proc.time()[["elapsed"]]
  for (scenario in names(scenarios)) {
    for (iteration in seq_len(simulations)) {
      position <- position + 1L
      iteration_seed <- as.integer(
        seed + match(scenario, names(scenarios)) * 100000L + iteration
      )
      set.seed(iteration_seed)
      block_effect <- rnorm(blocks, sd = 0.7)
      cell <- interaction(template$treatment, template$block, drop = TRUE)
      interaction_effect <- rnorm(nlevels(cell), sd = 1)
      response <- block_effect[template$block] + interaction_effect[cell] +
        scenarios[[scenario]] * (template$treatment == "A") +
        rnorm(nrow(template), sd = 0.7)
      design <- compile_permanova_design(
        ~ treatment + (1 | block) + (1 | treatment:block),
        template, stats::dist(response)
      )
      fit <- suppressWarnings(permanova_mixed_contrasts(
        design, contrast, permutations = permutations,
        seed = as.integer(iteration_seed + 50000000L)
      ))
      means <- stats::aggregate(
        response, list(treatment = template$treatment, block = template$block),
        mean
      )
      reference <- .scalar_mixed_contrast_reference(
        means$x, means$treatment, contrast, means$block
      )
      rows[[position]] <- data.frame(
        scenario = scenario, iteration = iteration,
        engine_f = fit$table$pseudo_f[[1L]],
        engine_p = fit$table$p_perm[[1L]],
        classical_f = reference$F,
        classical_p = stats::pf(reference$F, 1, reference$df,
                                lower.tail = FALSE),
        stringsAsFactors = FALSE
      )
    }
  }
  raw <- do.call(rbind, rows)
  summary <- do.call(rbind, lapply(names(scenarios), function(scenario) {
    values <- raw[raw$scenario == scenario, , drop = FALSE]
    do.call(rbind, lapply(c("engine", "classical"), function(method) {
      p <- values[[paste0(method, "_p")]]
      uncertainty <- .binomial_summary(sum(p <= alpha), length(p))
      data.frame(
        scenario = scenario, method = method, simulations = length(p),
        rejection_rate = uncertainty[["rejection_rate"]],
        monte_carlo_se = uncertainty[["monte_carlo_se"]],
        confidence_low = uncertainty[["confidence_low"]],
        confidence_high = uncertainty[["confidence_high"]],
        stringsAsFactors = FALSE
      )
    }))
  }))
  rownames(summary) <- NULL
  out <- list(
    summary = summary,
    max_engine_classical_f_difference = max(abs(
      raw$engine_f - raw$classical_f
    )),
    results = if (keep_results) raw else NULL,
    settings = list(
      simulations = simulations, permutations = as.integer(permutations),
      blocks = blocks, replicates = replicates, effect = effect,
      alpha = alpha, seed = seed
    ),
    elapsed_seconds = proc.time()[["elapsed"]] - start,
    interpretation = paste(
      "The engine p-value is empirical; the classical p-value is a Gaussian",
      "reference. The invariant is equality of the observed contrast F."
    )
  )
  class(out) <- "permanova_mixed_contrast_validation"
  out
}

#' @export
print.permanova_mixed_contrast_validation <- function(x, ...) {
  cat("Mixed-design contrast validation\n")
  print(x$summary, row.names = FALSE)
  cat("Maximum engine/classical F difference:",
      x$max_engine_classical_f_difference, "\n")
  cat(x$interpretation, "\n")
  invisible(x)
}

#' Validate global-null max-F multiplicity control
#'
#' Simulates an unequal-size, three-group one-factor design under the complete
#' equal-centroid null and under a location shift. It reports family rejection
#' for unadjusted pairwise tests, Holm adjustment, and the shared-permutation
#' max-F adjustment.
#'
#' @param simulations Number of simulations per scenario.
#' @param permutations Permutations per fit.
#' @param effect Location shift for group C in the power scenario.
#' @param alpha Rejection threshold.
#' @param seed Reproducibility seed.
#' @param keep_results Retain simulation-level results.
#' @return A global-null max-F validation report.
#' @export
simulate_maxT_validation <- function(
    simulations = 500L, permutations = 99L, effect = 0.8,
    alpha = 0.05, seed = 20260719L, keep_results = FALSE) {
  simulations <- .positive_integer(simulations, "simulations")
  if (length(permutations) != 1L || !is.finite(permutations) ||
      permutations < 0 || permutations != as.integer(permutations)) {
    stop("`permutations` must be one non-negative integer.", call. = FALSE)
  }
  if (!is.finite(effect) || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`effect` must be finite and `alpha` must lie between zero and one.",
         call. = FALSE)
  }
  group <- factor(rep(c("A", "B", "C"), c(4, 5, 6)),
                  levels = c("A", "B", "C"))
  data <- data.frame(group = group)
  scenarios <- c(null = 0, location_shift = effect)
  rows <- vector("list", length(scenarios) * simulations)
  position <- 0L
  start <- proc.time()[["elapsed"]]
  for (scenario in names(scenarios)) {
    for (iteration in seq_len(simulations)) {
      position <- position + 1L
      iteration_seed <- as.integer(
        seed + match(scenario, names(scenarios)) * 100000L + iteration
      )
      set.seed(iteration_seed)
      response <- rnorm(length(group)) +
        scenarios[[scenario]] * (group == "C")
      design <- compile_permanova_design(
        ~ group, data, stats::dist(response)
      )
      fit <- permanova_maxT_pairwise(
        design, permutations = permutations,
        seed = as.integer(iteration_seed + 50000000L)
      )
      raw <- fit$table$p_perm
      rows[[position]] <- data.frame(
        scenario = scenario, iteration = iteration,
        unadjusted_family_reject = any(raw <= alpha),
        holm_family_reject = any(stats::p.adjust(raw, "holm") <= alpha),
        maxT_family_reject = any(fit$table$p_maxT <= alpha),
        stringsAsFactors = FALSE
      )
    }
  }
  raw <- do.call(rbind, rows)
  methods <- c("unadjusted", "holm", "maxT")
  summary <- do.call(rbind, lapply(names(scenarios), function(scenario) {
    values <- raw[raw$scenario == scenario, , drop = FALSE]
    do.call(rbind, lapply(methods, function(method) {
      rejected <- values[[paste0(method, "_family_reject")]]
      uncertainty <- .binomial_summary(sum(rejected), length(rejected))
      data.frame(
        scenario = scenario, method = method,
        simulations = length(rejected),
        family_rejection_rate = uncertainty[["rejection_rate"]],
        monte_carlo_se = uncertainty[["monte_carlo_se"]],
        confidence_low = uncertainty[["confidence_low"]],
        confidence_high = uncertainty[["confidence_high"]],
        stringsAsFactors = FALSE
      )
    }))
  }))
  rownames(summary) <- NULL
  out <- list(
    summary = summary, results = if (keep_results) raw else NULL,
    settings = list(
      simulations = simulations, permutations = as.integer(permutations),
      group_sizes = c(A = 4L, B = 5L, C = 6L), effect = effect,
      alpha = alpha, seed = seed
    ),
    elapsed_seconds = proc.time()[["elapsed"]] - start,
    interpretation = paste(
      "The max-F rate is a complete-null weak-FWER diagnostic.",
      "It does not establish strong subset-null control."
    )
  )
  class(out) <- "permanova_maxT_validation"
  out
}

#' @export
print.permanova_maxT_validation <- function(x, ...) {
  cat("Global-null max-F validation\n")
  print(x$summary, row.names = FALSE)
  cat(x$interpretation, "\n")
  invisible(x)
}

#' Validate unequal-replication mixed planned contrasts
#'
#' Simulates a three-treatment random-intercept design in which every exact
#' cluster-size stratum contains all treatments. The A-minus-B contrast is
#' compared with a scalar GLH on equally weighted cluster centroids and with a
#' direct reconstruction of the restricted Freedman--Lane distribution.
#'
#' @param simulations Number of simulations per scenario.
#' @param permutations Permutations per fit.
#' @param effect A-minus-B location shift in the power scenario.
#' @param alpha Rejection threshold.
#' @param seed Reproducibility seed.
#' @param keep_results Retain simulation-level results.
#' @return An unbalanced mixed-contrast validation report.
#' @export
simulate_unbalanced_mixed_contrast_validation <- function(
    simulations = 500L, permutations = 99L, effect = 0.8,
    alpha = 0.05, seed = 20260720L, keep_results = FALSE) {
  simulations <- .positive_integer(simulations, "simulations")
  if (length(permutations) != 1L || !is.finite(permutations) ||
      permutations < 0 || permutations != as.integer(permutations)) {
    stop("`permutations` must be one non-negative integer.", call. = FALSE)
  }
  if (!is.finite(effect) || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`effect` must be finite and `alpha` must lie between zero and one.",
         call. = FALSE)
  }
  sizes <- rep(2:4, each = 3L)
  treatment_by_site <- factor(
    rep(c("A", "B", "C"), 3L), levels = c("A", "B", "C")
  )
  data <- data.frame(
    site = factor(rep(seq_along(sizes), sizes)),
    treatment = factor(rep(treatment_by_site, sizes),
                       levels = levels(treatment_by_site))
  )
  contrast <- c(A = 1, B = -1, C = 0)
  scenarios <- c(null = 0, location_shift = effect)
  rows <- vector("list", length(scenarios) * simulations)
  position <- 0L
  start <- proc.time()[["elapsed"]]
  for (scenario in names(scenarios)) {
    for (iteration in seq_len(simulations)) {
      position <- position + 1L
      data_seed <- as.integer(
        seed + match(scenario, names(scenarios)) * 100000L + iteration
      )
      permutation_seed <- as.integer(data_seed + 50000000L)
      set.seed(data_seed)
      site_effect <- rnorm(length(sizes), sd = 1)
      response <- site_effect[data$site] +
        scenarios[[scenario]] * (data$treatment == "A") +
        rnorm(nrow(data), sd = 0.7)
      design <- compile_permanova_design(
        ~ treatment + (1 | site), data, stats::dist(response)
      )
      fit <- suppressWarnings(permanova_mixed_contrasts(
        design, contrast, permutations = permutations,
        seed = permutation_seed, keep_permutations = TRUE
      ))
      structure <- .mixed_contrast_structure(design)
      units <- .prepare_mixed_contrast_units(structure)
      set.seed(permutation_seed)
      plan <- .mixed_contrast_plan(structure, permutations, 100000L)
      coefficient_row <- as.numeric(contrast %*% units$mean_rows)
      projector <- .linear_contrast_projector(
        units$X, coefficient_row, units$H_full, 1e-10
      )
      centroids <- vapply(
        structure$aggregated$blocks, function(index) mean(response[index]),
        numeric(1)
      )
      residual <- diag(length(centroids)) - projector$H_reduced
      direct_statistics <- apply(plan$orders, 1L, function(order) {
        permuted <- (projector$H_reduced +
          residual[order, , drop = FALSE]) %*% centroids
        sum(permuted * (projector$H %*% permuted)) /
          (sum(permuted * (units$H_denominator %*% permuted)) /
             units$denominator_df)
      })
      observed <- fit$table$pseudo_f[[1L]]
      direct_p <- if (plan$exact) {
        mean(direct_statistics >= observed - 1e-10)
      } else {
        (sum(direct_statistics >= observed - 1e-10) + 1) /
          (length(direct_statistics) + 1)
      }
      reference <- .scalar_mixed_contrast_reference(
        centroids, structure$treatment, contrast
      )
      rows[[position]] <- data.frame(
        scenario = scenario, iteration = iteration,
        engine_f = observed, engine_p = fit$table$p_perm[[1L]],
        direct_p = direct_p, classical_f = reference$F,
        classical_p = stats::pf(reference$F, 1, reference$df,
                                lower.tail = FALSE),
        permutation_difference = max(abs(
          fit$tests[[1L]]$permutation_statistics - direct_statistics
        )), stringsAsFactors = FALSE
      )
    }
  }
  raw <- do.call(rbind, rows)
  methods <- c("engine", "direct", "classical")
  summary <- do.call(rbind, lapply(names(scenarios), function(scenario) {
    values <- raw[raw$scenario == scenario, , drop = FALSE]
    do.call(rbind, lapply(methods, function(method) {
      p <- values[[paste0(method, "_p")]]
      uncertainty <- .binomial_summary(sum(p <= alpha), length(p))
      data.frame(
        scenario = scenario, method = method,
        simulations = length(p),
        rejection_rate = uncertainty[["rejection_rate"]],
        monte_carlo_se = uncertainty[["monte_carlo_se"]],
        confidence_low = uncertainty[["confidence_low"]],
        confidence_high = uncertainty[["confidence_high"]],
        stringsAsFactors = FALSE
      )
    }))
  }))
  rownames(summary) <- NULL
  out <- list(
    summary = summary,
    max_engine_classical_f_difference = max(abs(
      raw$engine_f - raw$classical_f
    )),
    max_engine_direct_p_difference = max(abs(
      raw$engine_p - raw$direct_p
    )),
    max_engine_direct_permutation_difference = max(
      raw$permutation_difference
    ),
    results = if (keep_results) raw else NULL,
    settings = list(
      simulations = simulations, permutations = as.integer(permutations),
      sizes = sizes, effect = effect, alpha = alpha, seed = seed,
      permutation_seed_offset = 50000000L
    ),
    elapsed_seconds = proc.time()[["elapsed"]] - start
  )
  class(out) <- "permanova_unbalanced_mixed_contrast_validation"
  out
}

#' @export
print.permanova_unbalanced_mixed_contrast_validation <- function(x, ...) {
  cat("Unbalanced mixed-contrast validation\n")
  print(x$summary, row.names = FALSE)
  cat("Maximum engine/classical F difference:",
      x$max_engine_classical_f_difference, "\n")
  cat("Maximum engine/direct permutation difference:",
      x$max_engine_direct_permutation_difference, "\n")
  invisible(x)
}

#' Validate explicitly covariate-adjusted contrasts
#'
#' Simulates an unbalanced fixed design with a categorical nuisance factor and
#' a continuous covariate. An adjusted A-minus-B contrast at an explicit
#' covariate value is compared with a scalar Gaussian GLH and with a direct
#' reconstruction of every sampled Freedman--Lane statistic.
#'
#' @param simulations Number of simulations per scenario.
#' @param permutations Permutations per fit.
#' @param effect A-minus-B location shift in the power scenario.
#' @param alpha Rejection threshold.
#' @param at Continuous-covariate reference value.
#' @param seed Reproducibility seed.
#' @param keep_results Retain simulation-level results.
#' @return A covariate-adjusted contrast validation report.
#' @export
simulate_adjusted_contrast_validation <- function(
    simulations = 500L, permutations = 99L, effect = 0.8,
    alpha = 0.05, at = 0.25, seed = 20260721L,
    keep_results = FALSE) {
  simulations <- .positive_integer(simulations, "simulations")
  if (length(permutations) != 1L || !is.finite(permutations) ||
      permutations < 0 || permutations != as.integer(permutations)) {
    stop("`permutations` must be one non-negative integer.", call. = FALSE)
  }
  if (!is.finite(effect) || !is.finite(alpha) || alpha <= 0 || alpha >= 1 ||
      length(at) != 1L || !is.finite(at)) {
    stop("`effect` and `at` must be finite and `alpha` must lie between zero and one.",
         call. = FALSE)
  }
  data <- expand.grid(
    replicate = seq_len(4L),
    treatment = factor(c("A", "B", "C"), levels = c("A", "B", "C")),
    habitat = factor(c("closed", "open"), levels = c("closed", "open"))
  )
  data <- data[-c(2L, 11L, 18L), , drop = FALSE]
  rownames(data) <- NULL
  data$gradient <- seq(-1.4, 1.6, length.out = nrow(data)) +
    c(A = -0.25, B = 0.05, C = 0.30)[data$treatment]
  contrast <- c(A = 1, B = -1, C = 0)
  scenarios <- c(null = 0, location_shift = effect)
  rows <- vector("list", length(scenarios) * simulations)
  position <- 0L
  start <- proc.time()[["elapsed"]]
  for (scenario in names(scenarios)) {
    for (iteration in seq_len(simulations)) {
      position <- position + 1L
      data_seed <- as.integer(
        seed + match(scenario, names(scenarios)) * 100000L + iteration
      )
      permutation_seed <- as.integer(data_seed + 50000000L)
      set.seed(data_seed)
      response <- 0.75 * data$gradient +
        0.35 * (data$habitat == "open") +
        scenarios[[scenario]] * (data$treatment == "A") +
        stats::rnorm(nrow(data))
      design <- compile_permanova_design(
        ~ treatment * gradient + habitat, data, stats::dist(response)
      )
      fit <- permanova_adjusted_contrasts(
        design, "treatment", contrast, at = c(gradient = at),
        permutations = permutations, seed = permutation_seed,
        keep_permutations = TRUE
      )
      means <- derive_adjusted_means(
        design, "treatment", at = c(gradient = at)
      )
      coefficient_row <- as.numeric(contrast %*% means$coefficient_rows)
      projector <- .linear_contrast_projector(
        design$X, coefficient_row, design$projections$fixed$H, 1e-10
      )
      set.seed(permutation_seed)
      plan <- .permutation_orders(nrow(data), permutations, 100000L)
      residual <- diag(nrow(data)) - projector$H_reduced
      denominator_H <- diag(nrow(data)) - design$projections$fixed$H
      denominator_df <- nrow(data) - design$estimability$fixed_rank
      direct_statistics <- apply(plan$orders, 1L, function(order) {
        permuted <- (projector$H_reduced +
          residual[order, , drop = FALSE]) %*% response
        sum(permuted * (projector$H %*% permuted)) /
          (sum(permuted * (denominator_H %*% permuted)) / denominator_df)
      })
      observed <- fit$table$pseudo_f[[1L]]
      direct_p <- if (plan$exact) {
        mean(direct_statistics >= observed - 1e-10)
      } else {
        (sum(direct_statistics >= observed - 1e-10) + 1) /
          (length(direct_statistics) + 1)
      }
      beta <- qr.solve(design$X, response)
      estimate <- sum(coefficient_row * beta)
      numerator <- estimate ^ 2 / as.numeric(
        coefficient_row %*% solve(crossprod(design$X), coefficient_row)
      )
      denominator <- sum((response - design$X %*% beta) ^ 2) /
        denominator_df
      classical_f <- numerator / denominator
      rows[[position]] <- data.frame(
        scenario = scenario, iteration = iteration,
        engine_f = observed, engine_p = fit$table$p_perm[[1L]],
        direct_p = direct_p, classical_f = classical_f,
        classical_p = stats::pf(classical_f, 1, denominator_df,
                                lower.tail = FALSE),
        permutation_difference = max(abs(
          fit$tests[[1L]]$permutation_statistics - direct_statistics
        )), stringsAsFactors = FALSE
      )
    }
  }
  raw <- do.call(rbind, rows)
  methods <- c("engine", "direct", "classical")
  summary <- do.call(rbind, lapply(names(scenarios), function(scenario) {
    values <- raw[raw$scenario == scenario, , drop = FALSE]
    do.call(rbind, lapply(methods, function(method) {
      p <- values[[paste0(method, "_p")]]
      uncertainty <- .binomial_summary(sum(p <= alpha), length(p))
      data.frame(
        scenario = scenario, method = method,
        simulations = length(p),
        rejection_rate = uncertainty[["rejection_rate"]],
        monte_carlo_se = uncertainty[["monte_carlo_se"]],
        confidence_low = uncertainty[["confidence_low"]],
        confidence_high = uncertainty[["confidence_high"]],
        stringsAsFactors = FALSE
      )
    }))
  }))
  rownames(summary) <- NULL
  out <- list(
    summary = summary,
    max_engine_classical_f_difference = max(abs(
      raw$engine_f - raw$classical_f
    )),
    max_engine_direct_p_difference = max(abs(
      raw$engine_p - raw$direct_p
    )),
    max_engine_direct_permutation_difference = max(
      raw$permutation_difference
    ),
    results = if (keep_results) raw else NULL,
    settings = list(
      simulations = simulations, permutations = as.integer(permutations),
      effect = effect, alpha = alpha, at = at, seed = seed,
      permutation_seed_offset = 50000000L
    ),
    elapsed_seconds = proc.time()[["elapsed"]] - start
  )
  class(out) <- "permanova_adjusted_contrast_validation"
  out
}

#' @export
print.permanova_adjusted_contrast_validation <- function(x, ...) {
  cat("Covariate-adjusted contrast validation\n")
  print(x$summary, row.names = FALSE)
  cat("Maximum engine/classical F difference:",
      x$max_engine_classical_f_difference, "\n")
  cat("Maximum engine/direct permutation difference:",
      x$max_engine_direct_permutation_difference, "\n")
  invisible(x)
}

#' Validate mixed-design global-null max-F inference
#'
#' Simulates an unequally replicated three-treatment random-intercept design
#' whose exact-size strata each contain every treatment. The complete shared
#' statistic matrix is reconstructed from scalar unit centroids, and family
#' rejection is reported for unadjusted, Holm, and weak-FWER max-F procedures.
#'
#' @param simulations Number of simulations per scenario.
#' @param permutations Permutations per fit.
#' @param effect Group-C location shift in the power scenario.
#' @param alpha Family rejection threshold.
#' @param seed Reproducibility seed.
#' @param keep_results Retain simulation-level results.
#' @return A mixed global-null max-F validation report.
#' @export
simulate_mixed_maxT_validation <- function(
    simulations = 500L, permutations = 99L, effect = 0.8,
    alpha = 0.05, seed = 20260722L, keep_results = FALSE) {
  simulations <- .positive_integer(simulations, "simulations")
  if (length(permutations) != 1L || !is.finite(permutations) ||
      permutations < 0 || permutations != as.integer(permutations)) {
    stop("`permutations` must be one non-negative integer.", call. = FALSE)
  }
  if (!is.finite(effect) || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`effect` must be finite and `alpha` must lie between zero and one.",
         call. = FALSE)
  }
  sizes <- rep(2:4, each = 3L)
  treatment_by_site <- factor(
    rep(c("A", "B", "C"), 3L), levels = c("A", "B", "C")
  )
  data <- data.frame(
    site = factor(rep(seq_along(sizes), sizes)),
    treatment = factor(rep(treatment_by_site, sizes),
                       levels = levels(treatment_by_site))
  )
  scenarios <- c(null = 0, location_shift = effect)
  rows <- vector("list", length(scenarios) * simulations)
  position <- 0L
  start <- proc.time()[["elapsed"]]
  for (scenario in names(scenarios)) {
    for (iteration in seq_len(simulations)) {
      position <- position + 1L
      data_seed <- as.integer(
        seed + match(scenario, names(scenarios)) * 100000L + iteration
      )
      permutation_seed <- as.integer(data_seed + 50000000L)
      set.seed(data_seed)
      site_effect <- stats::rnorm(length(sizes), sd = 1)
      response <- site_effect[data$site] +
        scenarios[[scenario]] * (data$treatment == "C") +
        stats::rnorm(nrow(data), sd = 0.7)
      design <- compile_permanova_design(
        ~ treatment + (1 | site), data, stats::dist(response)
      )
      fit <- permanova_mixed_maxT_pairwise(
        design, permutations = permutations, seed = permutation_seed,
        keep_permutations = TRUE
      )
      structure <- .mixed_contrast_structure(design)
      units <- .prepare_mixed_contrast_units(structure)
      set.seed(permutation_seed)
      plan <- .mixed_contrast_plan(structure, permutations, 100000L)
      projectors <- lapply(
        seq_len(nrow(fit$contrast_matrix)), function(index) {
          coefficient_row <- as.numeric(
            fit$contrast_matrix[index, ] %*% units$mean_rows
          )
          .linear_contrast_projector(
            units$X, coefficient_row, units$H_full, 1e-10
          )$H
        }
      )
      centroids <- vapply(
        structure$aggregated$blocks, function(index) mean(response[index]),
        numeric(1)
      )
      scalar_statistics <- t(apply(plan$orders, 1L, function(order) {
        permuted <- centroids[order]
        denominator <- sum(
          permuted * (units$H_denominator %*% permuted)
        ) / units$denominator_df
        vapply(projectors, function(H) {
          sum(permuted * (H %*% permuted)) / denominator
        }, numeric(1))
      }))
      scalar_observed <- vapply(projectors, function(H) {
        numerator <- sum(centroids * (H %*% centroids))
        denominator <- sum(
          centroids * (units$H_denominator %*% centroids)
        ) / units$denominator_df
        numerator / denominator
      }, numeric(1))
      maxima <- apply(scalar_statistics, 1L, max)
      direct_raw <- direct_max <- numeric(length(scalar_observed))
      for (index in seq_along(scalar_observed)) {
        raw_extreme <- sum(
          scalar_statistics[, index] >= scalar_observed[[index]] - 1e-10
        )
        max_extreme <- sum(maxima >= scalar_observed[[index]] - 1e-10)
        if (plan$exact) {
          direct_raw[[index]] <- raw_extreme / nrow(scalar_statistics)
          direct_max[[index]] <- max_extreme / nrow(scalar_statistics)
        } else {
          direct_raw[[index]] <- (raw_extreme + 1) /
            (nrow(scalar_statistics) + 1)
          direct_max[[index]] <- (max_extreme + 1) /
            (nrow(scalar_statistics) + 1)
        }
      }
      raw <- fit$table$p_perm
      rows[[position]] <- data.frame(
        scenario = scenario, iteration = iteration,
        unadjusted_family_reject = any(raw <= alpha),
        holm_family_reject = any(stats::p.adjust(raw, "holm") <= alpha),
        maxT_family_reject = any(fit$table$p_maxT <= alpha),
        f_difference = max(abs(fit$table$pseudo_f - scalar_observed)),
        statistic_difference = max(abs(
          fit$permutation$statistic_matrix - scalar_statistics
        )),
        p_difference = max(abs(c(
          fit$table$p_perm - direct_raw,
          fit$table$p_maxT - direct_max
        ))), stringsAsFactors = FALSE
      )
    }
  }
  raw <- do.call(rbind, rows)
  methods <- c("unadjusted", "holm", "maxT")
  summary <- do.call(rbind, lapply(names(scenarios), function(scenario) {
    values <- raw[raw$scenario == scenario, , drop = FALSE]
    do.call(rbind, lapply(methods, function(method) {
      rejected <- values[[paste0(method, "_family_reject")]]
      uncertainty <- .binomial_summary(sum(rejected), length(rejected))
      data.frame(
        scenario = scenario, method = method,
        simulations = length(rejected),
        family_rejection_rate = uncertainty[["rejection_rate"]],
        monte_carlo_se = uncertainty[["monte_carlo_se"]],
        confidence_low = uncertainty[["confidence_low"]],
        confidence_high = uncertainty[["confidence_high"]],
        stringsAsFactors = FALSE
      )
    }))
  }))
  rownames(summary) <- NULL
  out <- list(
    summary = summary,
    max_engine_scalar_f_difference = max(raw$f_difference),
    max_engine_scalar_statistic_difference = max(raw$statistic_difference),
    max_engine_scalar_p_difference = max(raw$p_difference),
    results = if (keep_results) raw else NULL,
    settings = list(
      simulations = simulations, permutations = as.integer(permutations),
      sizes = sizes, effect = effect, alpha = alpha, seed = seed,
      permutation_seed_offset = 50000000L
    ),
    elapsed_seconds = proc.time()[["elapsed"]] - start,
    interpretation = paste(
      "The max-F rate is a complete-treatment-null weak-FWER diagnostic.",
      "It does not establish strong subset-null control."
    )
  )
  class(out) <- "permanova_mixed_maxT_validation"
  out
}

#' @export
print.permanova_mixed_maxT_validation <- function(x, ...) {
  cat("Mixed-design global-null max-F validation\n")
  print(x$summary, row.names = FALSE)
  cat("Maximum engine/scalar statistic difference:",
      x$max_engine_scalar_statistic_difference, "\n")
  cat(x$interpretation, "\n")
  invisible(x)
}

#' Validate lower-level unbalanced nested inference
#'
#' Simulates child-level and observation-level treatments in one unbalanced
#' parent/child hierarchy. Each distance result is compared with a scalar unit-
#' level general linear model and with direct reconstruction of the restricted
#' Freedman--Lane permutation distribution.
#'
#' @param simulations Number of simulations per structure and scenario.
#' @param permutations Permutations per fit.
#' @param effect Treatment location shift in the power scenario.
#' @param alpha Rejection threshold.
#' @param seed Reproducibility seed.
#' @param keep_results Retain simulation-level results.
#' @return A lower-level unbalanced nested validation report.
#' @export
simulate_unbalanced_nested_lower_validation <- function(
    simulations = 500L, permutations = 99L, effect = 0.8,
    alpha = 0.05, seed = 20260723L, keep_results = FALSE) {
  simulations <- .positive_integer(simulations, "simulations")
  if (length(permutations) != 1L || !is.finite(permutations) ||
      permutations < 0 || permutations != as.integer(permutations)) {
    stop("`permutations` must be one non-negative integer.", call. = FALSE)
  }
  if (!is.finite(effect) || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`effect` must be finite and `alpha` must lie between zero and one.",
         call. = FALSE)
  }
  child_rows <- list()
  child_id <- 0L
  for (parent in 1:4) {
    for (size in 2:3) {
      for (treatment in c("A", "B")) {
        child_id <- child_id + 1L
        child_rows[[length(child_rows) + 1L]] <- data.frame(
          parent = parent, child = child_id, treatment = treatment,
          observation = seq_len(size)
        )
      }
    }
  }
  child_data <- do.call(rbind, child_rows)
  observation_parent <- c(1, 1, 2, 2, 2, 3, 3, 3, 3)
  observation_data <- do.call(rbind, lapply(
    seq_along(observation_parent), function(child) {
      data.frame(
        parent = observation_parent[[child]], child = child,
        treatment = rep(c("A", "B"), 2L), observation = 1:4
      )
    }
  ))
  prepare <- function(data) {
    data$parent <- factor(data$parent)
    data$child <- factor(data$child)
    data$treatment <- factor(data$treatment, levels = c("A", "B"))
    data
  }
  templates <- list(
    child = prepare(child_data), observation = prepare(observation_data)
  )
  scenarios <- c(null = 0, location_shift = effect)
  rows <- vector(
    "list", length(templates) * length(scenarios) * simulations
  )
  position <- 0L
  start <- proc.time()[["elapsed"]]
  for (structure_name in names(templates)) {
    data <- templates[[structure_name]]
    for (scenario in names(scenarios)) {
      for (iteration in seq_len(simulations)) {
        position <- position + 1L
        data_seed <- as.integer(
          seed + match(structure_name, names(templates)) * 1000000L +
            match(scenario, names(scenarios)) * 100000L + iteration
        )
        permutation_seed <- as.integer(data_seed + 50000000L)
        set.seed(data_seed)
        parent_effect <- stats::rnorm(nlevels(data$parent), sd = 0.8)
        child_effect <- stats::rnorm(nlevels(data$child), sd = 0.7)
        response <- parent_effect[data$parent] + child_effect[data$child] +
          scenarios[[scenario]] * (data$treatment == "B") +
          stats::rnorm(nrow(data), sd = 0.6)
        design <- compile_permanova_design(
          ~ treatment + (1 | parent / child), data, stats::dist(response)
        )
        fit <- permanova_unbalanced_nested_lower(
          design, level = structure_name, permutations = permutations,
          seed = permutation_seed, keep_permutations = TRUE
        )
        validated <- .validate_unbalanced_nested_lower(
          design, level = structure_name
        )
        set.seed(permutation_seed)
        plan <- .mixed_contrast_plan(
          validated$structure, permutations, 100000L
        )
        unit_response <- vapply(
          validated$structure$aggregated$blocks,
          function(index) mean(response[index]), numeric(1)
        )
        residual <- diag(length(unit_response)) - validated$H_reduced
        direct_statistics <- apply(plan$orders, 1L, function(order) {
          permuted <- (validated$H_reduced +
            residual[order, , drop = FALSE]) %*% unit_response
          (sum(permuted * (validated$H_numerator %*% permuted)) /
             validated$numerator_df) /
            (sum(permuted *
                   (validated$units$H_denominator %*% permuted)) /
               validated$units$denominator_df)
        })
        observed <- fit$table$pseudo_f[[1L]]
        direct_p <- if (plan$exact) {
          mean(direct_statistics >= observed - 1e-10)
        } else {
          (sum(direct_statistics >= observed - 1e-10) + 1) /
            (length(direct_statistics) + 1)
        }
        scalar_data <- validated$units$data
        scalar_data$response <- unit_response
        scalar_fit <- stats::lm(
          response ~ .stratum + .treatment, data = scalar_data
        )
        scalar_table <- stats::anova(scalar_fit)
        classical_f <- scalar_table[".treatment", "F value"]
        classical_p <- scalar_table[".treatment", "Pr(>F)"]
        rows[[position]] <- data.frame(
          structure = structure_name, scenario = scenario,
          iteration = iteration, engine_f = observed,
          engine_p = fit$table$p_perm[[1L]], direct_p = direct_p,
          classical_f = classical_f, classical_p = classical_p,
          permutation_difference = max(abs(
            fit$test$permutation_statistics - direct_statistics
          )), stringsAsFactors = FALSE
        )
      }
    }
  }
  raw <- do.call(rbind, rows)
  methods <- c("engine", "direct", "classical")
  summary <- do.call(rbind, lapply(names(templates), function(structure_name) {
    do.call(rbind, lapply(names(scenarios), function(scenario) {
      values <- raw[
        raw$structure == structure_name & raw$scenario == scenario,
        , drop = FALSE
      ]
      do.call(rbind, lapply(methods, function(method) {
        p <- values[[paste0(method, "_p")]]
        uncertainty <- .binomial_summary(sum(p <= alpha), length(p))
        data.frame(
          structure = structure_name, scenario = scenario, method = method,
          simulations = length(p),
          rejection_rate = uncertainty[["rejection_rate"]],
          monte_carlo_se = uncertainty[["monte_carlo_se"]],
          confidence_low = uncertainty[["confidence_low"]],
          confidence_high = uncertainty[["confidence_high"]],
          stringsAsFactors = FALSE
        )
      }))
    }))
  }))
  rownames(summary) <- NULL
  out <- list(
    summary = summary,
    max_engine_classical_f_difference = max(abs(
      raw$engine_f - raw$classical_f
    )),
    max_engine_direct_p_difference = max(abs(
      raw$engine_p - raw$direct_p
    )),
    max_engine_direct_permutation_difference = max(
      raw$permutation_difference
    ),
    results = if (keep_results) raw else NULL,
    settings = list(
      simulations = simulations, permutations = as.integer(permutations),
      effect = effect, alpha = alpha, seed = seed,
      permutation_seed_offset = 50000000L
    ),
    elapsed_seconds = proc.time()[["elapsed"]] - start
  )
  out$agreement <- do.call(rbind, lapply(names(templates), function(name) {
    values <- raw[raw$structure == name, , drop = FALSE]
    data.frame(
      structure = name,
      max_f_difference = max(abs(values$engine_f - values$classical_f)),
      max_p_difference = max(abs(values$engine_p - values$direct_p)),
      max_permutation_difference = max(values$permutation_difference),
      stringsAsFactors = FALSE
    )
  }))
  rownames(out$agreement) <- NULL
  class(out) <- "permanova_unbalanced_nested_lower_validation"
  out
}

#' @export
print.permanova_unbalanced_nested_lower_validation <- function(x, ...) {
  cat("Lower-level unbalanced nested validation\n")
  print(x$summary, row.names = FALSE)
  cat("Maximum engine/classical F difference:",
      x$max_engine_classical_f_difference, "\n")
  cat("Maximum engine/direct permutation difference:",
      x$max_engine_direct_permutation_difference, "\n")
  invisible(x)
}

#' Validate unbalanced crossed-nested hybrid inference
#'
#' Simulates one treatment crossed with an unequal region/site hierarchy with
#' one missing treatment-by-region profile and unequal replicated site cells.
#' The distance engine is compared with the scalar equal-site-weight profile
#' model and direct reconstruction of every sampled profile permutation.
#'
#' @param simulations Number of simulations per scenario.
#' @param permutations Permutations per fit.
#' @param effect Treatment location shift in the power scenario.
#' @param alpha Rejection threshold.
#' @param seed Reproducibility seed.
#' @param keep_results Retain simulation-level results.
#' @return An unbalanced hybrid validation report.
#' @export
simulate_unbalanced_hybrid_validation <- function(
    simulations = 500L, permutations = 99L, effect = 0.8,
    alpha = 0.05, seed = 20260724L, keep_results = FALSE) {
  simulations <- .positive_integer(simulations, "simulations")
  if (length(permutations) != 1L || !is.finite(permutations) ||
      permutations < 0 || permutations != as.integer(permutations)) {
    stop("`permutations` must be one non-negative integer.", call. = FALSE)
  }
  if (!is.finite(effect) || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`effect` must be finite and `alpha` must lie between zero and one.",
         call. = FALSE)
  }
  sites_per_region <- c(R1 = 2L, R2 = 3L, R3 = 2L, R4 = 3L)
  treatments <- c("A", "B", "C")
  rows <- list()
  position <- 0L
  for (region in names(sites_per_region)) {
    region_index <- match(region, names(sites_per_region))
    for (treatment in treatments) {
      if (region == "R4" && treatment == "C") next
      for (site_index in seq_len(sites_per_region[[region]])) {
        cell_size <- 1L + ((region_index + site_index) %% 2L)
        position <- position + 1L
        rows[[position]] <- data.frame(
          treatment = treatment, region = region,
          site = paste(region, paste0("S", site_index), sep = "_"),
          replicate = seq_len(cell_size), stringsAsFactors = FALSE
        )
      }
    }
  }
  template <- do.call(rbind, rows)
  template$treatment <- factor(template$treatment, levels = treatments)
  template$region <- factor(template$region, levels = names(sites_per_region))
  template$site <- factor(template$site)
  region_site <- interaction(template$region, template$site, drop = TRUE)
  treatment_region <- interaction(
    template$treatment, template$region, drop = TRUE
  )
  treatment_site <- interaction(
    template$treatment, template$region, template$site, drop = TRUE
  )
  scenarios <- c(null = 0, location_shift = effect)
  results <- vector("list", length(scenarios) * simulations)
  result_index <- 0L
  start <- proc.time()[["elapsed"]]
  for (scenario in names(scenarios)) {
    for (iteration in seq_len(simulations)) {
      result_index <- result_index + 1L
      data_seed <- as.integer(
        seed + match(scenario, names(scenarios)) * 1000000L + iteration
      )
      permutation_seed <- as.integer(data_seed + 50000000L)
      set.seed(data_seed)
      response <-
        stats::rnorm(nlevels(template$region), sd = 0.8)[template$region] +
        stats::rnorm(nlevels(region_site), sd = 0.6)[region_site] +
        stats::rnorm(nlevels(treatment_region), sd = 0.7)[treatment_region] +
        stats::rnorm(nlevels(treatment_site), sd = 0.4)[treatment_site] +
        scenarios[[scenario]] * (template$treatment == "C") +
        stats::rnorm(nrow(template), sd = 0.5)
      design <- compile_permanova_design(
        ~ treatment + (1 | region/site) +
          (1 | treatment:region) + (1 | treatment:region:site),
        template, stats::dist(response)
      )
      fit <- permanova_unbalanced_hybrid(
        design, permutations = permutations, seed = permutation_seed,
        keep_permutations = TRUE
      )
      structure <- .validate_unbalanced_hybrid(design)
      profile_response <- drop(
        structure$aggregated$membership %*% response
      )
      scalar_frame <- data.frame(
        response = profile_response, region = structure$region,
        treatment = structure$treatment
      )
      scalar_fit <- stats::lm(response ~ region + treatment, scalar_frame)
      scalar_table <- stats::anova(scalar_fit)
      classical_f <- scalar_table["treatment", "F value"]
      classical_p <- scalar_table["treatment", "Pr(>F)"]
      set.seed(permutation_seed)
      plan <- .stratified_permutation_plan(
        structure$exchangeability_strata, permutations, 100000L
      )
      direct_statistics <- apply(plan$orders, 1L, function(order) {
        permuted <- profile_response[order]
        (sum(permuted * (structure$H_numerator %*% permuted)) /
           structure$numerator_df) /
          (sum(permuted * (structure$H_denominator %*% permuted)) /
             structure$denominator_df)
      })
      observed <- fit$test$statistic
      direct_p <- if (plan$exact) {
        mean(direct_statistics >= observed - 1e-10)
      } else {
        (sum(direct_statistics >= observed - 1e-10) + 1) /
          (length(direct_statistics) + 1)
      }
      results[[result_index]] <- data.frame(
        scenario = scenario, iteration = iteration,
        engine_f = observed, engine_p = fit$test$p_value,
        direct_p = direct_p, classical_f = classical_f,
        classical_p = classical_p,
        permutation_difference = max(abs(
          fit$test$permutation_statistics - direct_statistics
        )), stringsAsFactors = FALSE
      )
    }
  }
  raw <- do.call(rbind, results)
  methods <- c("engine", "direct", "classical")
  summary <- do.call(rbind, lapply(names(scenarios), function(scenario) {
    values <- raw[raw$scenario == scenario, , drop = FALSE]
    do.call(rbind, lapply(methods, function(method) {
      p <- values[[paste0(method, "_p")]]
      uncertainty <- .binomial_summary(sum(p <= alpha), length(p))
      data.frame(
        scenario = scenario, method = method, simulations = length(p),
        rejection_rate = uncertainty[["rejection_rate"]],
        monte_carlo_se = uncertainty[["monte_carlo_se"]],
        confidence_low = uncertainty[["confidence_low"]],
        confidence_high = uncertainty[["confidence_high"]],
        stringsAsFactors = FALSE
      )
    }))
  }))
  rownames(summary) <- NULL
  out <- list(
    summary = summary,
    max_engine_classical_f_difference = max(abs(
      raw$engine_f - raw$classical_f
    )),
    max_engine_direct_p_difference = max(abs(
      raw$engine_p - raw$direct_p
    )),
    max_engine_direct_permutation_difference = max(
      raw$permutation_difference
    ),
    results = if (keep_results) raw else NULL,
    settings = list(
      simulations = simulations, permutations = as.integer(permutations),
      effect = effect, alpha = alpha, seed = seed,
      permutation_seed_offset = 50000000L
    ),
    elapsed_seconds = proc.time()[["elapsed"]] - start
  )
  class(out) <- "permanova_unbalanced_hybrid_validation"
  out
}

#' @export
print.permanova_unbalanced_hybrid_validation <- function(x, ...) {
  cat("Unbalanced crossed-nested hybrid validation\n")
  print(x$summary, row.names = FALSE)
  cat("Maximum engine/classical F difference:",
      x$max_engine_classical_f_difference, "\n")
  cat("Maximum engine/direct permutation difference:",
      x$max_engine_direct_permutation_difference, "\n")
  invisible(x)
}

#' Validate balanced two-condition random-slope inference
#'
#' Compares the distance engine with paired-difference scalar F tests, direct
#' whole-subject sign-flip reconstruction, and scalar method-of-moments random
#' intercept--slope covariance estimates. A centred asymmetric null is included
#' as a labelled sensitivity outside the central-symmetry guarantee.
#'
#' @param simulations Number of simulations per scenario.
#' @param permutations Sign flips per fit.
#' @param effect Fixed mean difference in the power scenario.
#' @param alpha Rejection threshold.
#' @param seed Reproducibility seed.
#' @param keep_results Retain simulation-level results.
#' @return A random-slope validation report.
#' @export
simulate_random_slope_validation <- function(
    simulations = 500L, permutations = 99L, effect = 0.8,
    alpha = 0.05, seed = 20260725L, keep_results = FALSE) {
  simulations <- .positive_integer(simulations, "simulations")
  if (length(permutations) != 1L || !is.finite(permutations) ||
      permutations < 0 || permutations != as.integer(permutations)) {
    stop("`permutations` must be one non-negative integer.", call. = FALSE)
  }
  if (!is.finite(effect) || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`effect` must be finite and `alpha` must lie between zero and one.",
         call. = FALSE)
  }
  subjects <- 8L
  replicates <- 2L
  data <- expand.grid(
    replicate = seq_len(replicates),
    condition = factor(c("before", "after"), levels = c("before", "after")),
    subject = factor(seq_len(subjects))
  )
  scenarios <- c(
    symmetric_null = 0, location_shift = effect,
    asymmetric_zero_mean = 0
  )
  rows <- vector("list", length(scenarios) * simulations)
  row_index <- 0L
  start <- proc.time()[["elapsed"]]
  for (scenario in names(scenarios)) {
    for (iteration in seq_len(simulations)) {
      row_index <- row_index + 1L
      data_seed <- as.integer(
        seed + match(scenario, names(scenarios)) * 1000000L + iteration
      )
      permutation_seed <- as.integer(data_seed + 50000000L)
      set.seed(data_seed)
      latent_intercept <- stats::rnorm(subjects)
      random_intercept <- 0.8 * latent_intercept
      random_slope <- if (scenario == "asymmetric_zero_mean") {
        0.7 * (stats::rexp(subjects) - 1)
      } else {
        0.7 * (
          0.45 * latent_intercept + sqrt(1 - 0.45^2) * stats::rnorm(subjects)
        )
      }
      response <- random_intercept[data$subject] +
        (random_slope[data$subject] + scenarios[[scenario]]) *
          (data$condition == "after") +
        stats::rnorm(nrow(data), sd = 0.45)
      design <- compile_permanova_design(
        ~ condition + (1 + condition | subject),
        data, stats::dist(response)
      )
      fit <- permanova_random_slope(
        design, permutations = permutations, seed = permutation_seed,
        keep_permutations = TRUE
      )
      structure <- .validate_random_slope_design(design)
      difference <- drop(structure$difference_membership %*% response)
      scalar_f <- mean(difference)^2 /
        (stats::var(difference) / length(difference))
      set.seed(permutation_seed)
      plan <- .sign_flip_plan(subjects, permutations, 100000L)
      direct_statistics <- apply(plan$signs, 1L, function(sign) {
        value <- sign * difference
        mean(value)^2 / (stats::var(value) / length(value))
      })
      direct_p <- if (plan$exact) {
        mean(direct_statistics >= fit$test$statistic - 1e-10)
      } else {
        (sum(direct_statistics >= fit$test$statistic - 1e-10) + 1) /
          (length(direct_statistics) + 1)
      }
      midpoint <- drop(structure$midpoint_membership %*% response)
      coefficients <- cbind(midpoint = midpoint, difference = difference)
      cell <- interaction(data$subject, data$condition)
      residual <- response - stats::ave(response, cell)
      residual_ms <- sum(residual^2) / (nrow(data) - nlevels(cell))
      covariance_reference <- stats::cov(coefficients)
      covariance_reference[1L, 1L] <-
        covariance_reference[1L, 1L] - residual_ms / (2 * replicates)
      covariance_reference[2L, 2L] <-
        covariance_reference[2L, 2L] - 2 * residual_ms / replicates
      dimnames(covariance_reference) <- dimnames(fit$covariance$unrestricted)
      rows[[row_index]] <- data.frame(
        scenario = scenario, iteration = iteration,
        engine_f = fit$test$statistic, scalar_f = scalar_f,
        engine_p = fit$test$p_value, direct_p = direct_p,
        parametric_p = fit$test$p_value_parametric,
        permutation_difference = max(abs(
          fit$test$permutation_statistics - direct_statistics
        )),
        covariance_difference = max(abs(
          fit$covariance$unrestricted - covariance_reference
        )), stringsAsFactors = FALSE
      )
    }
  }
  raw <- do.call(rbind, rows)
  methods <- c("engine", "direct", "parametric")
  summary <- do.call(rbind, lapply(names(scenarios), function(scenario) {
    values <- raw[raw$scenario == scenario, , drop = FALSE]
    do.call(rbind, lapply(methods, function(method) {
      p <- values[[paste0(method, "_p")]]
      uncertainty <- .binomial_summary(sum(p <= alpha), length(p))
      data.frame(
        scenario = scenario, method = method, simulations = length(p),
        rejection_rate = uncertainty[["rejection_rate"]],
        monte_carlo_se = uncertainty[["monte_carlo_se"]],
        confidence_low = uncertainty[["confidence_low"]],
        confidence_high = uncertainty[["confidence_high"]],
        guarantee = if (scenario == "asymmetric_zero_mean") {
          "sensitivity outside central symmetry"
        } else if (method == "parametric") {
          "Gaussian paired-difference reference"
        } else "central-symmetry sign-flip inference",
        stringsAsFactors = FALSE
      )
    }))
  }))
  rownames(summary) <- NULL
  out <- list(
    summary = summary,
    max_engine_scalar_f_difference = max(abs(raw$engine_f - raw$scalar_f)),
    max_engine_direct_p_difference = max(abs(raw$engine_p - raw$direct_p)),
    max_engine_direct_permutation_difference = max(
      raw$permutation_difference
    ),
    max_engine_scalar_covariance_difference = max(
      raw$covariance_difference
    ),
    results = if (keep_results) raw else NULL,
    settings = list(
      simulations = simulations, permutations = as.integer(permutations),
      effect = effect, alpha = alpha, seed = seed,
      subjects = subjects, replicates = replicates,
      permutation_seed_offset = 50000000L
    ),
    interpretation = paste(
      "Calibration claims apply only to centrally symmetric subject",
      "difference vectors; the asymmetric scenario is diagnostic only."
    ),
    elapsed_seconds = proc.time()[["elapsed"]] - start
  )
  class(out) <- "permanova_random_slope_validation"
  out
}

#' @export
print.permanova_random_slope_validation <- function(x, ...) {
  cat("Balanced two-condition random-slope validation\n")
  print(x$summary, row.names = FALSE)
  cat("Maximum engine/scalar F difference:",
      x$max_engine_scalar_f_difference, "\n")
  cat("Maximum engine/direct sign-flip difference:",
      x$max_engine_direct_permutation_difference, "\n")
  cat("Interpretation:", x$interpretation, "\n")
  invisible(x)
}

#' Validate strong closed pairwise permutation inference
#'
#' Audits complete- and partial-null familywise rejection, reconstructs every
#' sampled intersection maximum directly from the retained label assignments,
#' and verifies that the partial `A = B` intersection never moves observations
#' from the non-null `C` group.
#'
#' @param simulations Number of simulations per scenario.
#' @param permutations Label assignments per intersection.
#' @param effect Location shift for the non-null group in the partial-null
#'   scenario.
#' @param alpha Familywise rejection threshold.
#' @param seed Reproducibility seed.
#' @param keep_results Retain simulation-level results.
#' @return A strong-pairwise validation report.
#' @export
simulate_strong_pairwise_validation <- function(
    simulations = 500L, permutations = 99L, effect = 1.2,
    alpha = 0.05, seed = 20260726L, keep_results = FALSE) {
  simulations <- .positive_integer(simulations, "simulations")
  if (length(permutations) != 1L || !is.finite(permutations) ||
      permutations < 1 || permutations != as.integer(permutations)) {
    stop("`permutations` must be one positive integer.", call. = FALSE)
  }
  if (!is.finite(effect) || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`effect` must be finite and `alpha` must lie between zero and one.",
         call. = FALSE)
  }
  per_group <- 6L
  data <- data.frame(
    group = factor(rep(c("A", "B", "C"), each = per_group))
  )
  scenarios <- c(complete_null = 0, partial_null = effect)
  rows <- vector("list", length(scenarios) * simulations)
  row_index <- 0L
  start <- proc.time()[["elapsed"]]
  for (scenario in names(scenarios)) {
    for (iteration in seq_len(simulations)) {
      row_index <- row_index + 1L
      data_seed <- as.integer(
        seed + match(scenario, names(scenarios)) * 1000000L + iteration
      )
      permutation_seed <- as.integer(data_seed + 50000000L)
      set.seed(data_seed)
      response <- stats::rnorm(nrow(data)) +
        scenarios[[scenario]] * (data$group == "C")
      design <- compile_permanova_design(
        ~ group, data, stats::dist(response)
      )
      fit <- permanova_pairwise_strong(
        design, permutations = permutations, seed = permutation_seed,
        keep_permutations = TRUE
      )
      structure <- .validate_contrast_design(design)
      projectors <- lapply(seq_len(nrow(fit$contrast_matrix)), function(index) {
        .contrast_projector(
          fit$contrast_matrix[index, ], structure$group, structure$counts,
          structure$H_full, 1e-10
        )$H
      })
      direct_intersection_p <- numeric(length(fit$intersection_details))
      statistic_difference <- numeric(length(fit$intersection_details))
      for (intersection in seq_along(fit$intersection_details)) {
        detail <- fit$intersection_details[[intersection]]
        direct_statistics <- vapply(seq_len(nrow(detail$orders)), function(row) {
          order <- detail$orders[row, ]
          G <- design$G[order, order, drop = FALSE]
          denominator <- trace_hgh(structure$H_denominator, G) /
            structure$denominator_df
          max(vapply(detail$subset, function(index) {
            trace_hgh(projectors[[index]], G) / denominator
          }, numeric(1)))
        }, numeric(1))
        statistic_difference[[intersection]] <- max(abs(
          direct_statistics - detail$statistics
        ))
        extreme <- sum(direct_statistics >= detail$statistic - 1e-10)
        direct_intersection_p[[intersection]] <- if (detail$exact) {
          extreme / length(direct_statistics)
        } else (extreme + 1) / (length(direct_statistics) + 1)
      }
      direct_closed <- vapply(seq_len(nrow(fit$table)), function(index) {
        containing <- vapply(fit$intersection_details, function(detail) {
          index %in% detail$subset
        }, logical(1))
        max(direct_intersection_p[containing])
      }, numeric(1))
      ab <- fit$intersection_details[[which(vapply(
        fit$intersection_details,
        function(detail) identical(detail$hypotheses, "A - B"), logical(1)
      ))]]
      c_positions <- which(data$group == "C")
      preservation_violations <- sum(
        ab$orders[, c_positions, drop = FALSE] !=
          matrix(c_positions, nrow(ab$orders), length(c_positions),
                 byrow = TRUE)
      )
      rows[[row_index]] <- data.frame(
        scenario = scenario, iteration = iteration,
        family_reject = any(fit$table$p_closed <= alpha),
        true_null_reject = if (scenario == "complete_null") {
          any(fit$table$p_closed <= alpha)
        } else fit$table$p_closed[fit$table$contrast == "A - B"] <= alpha,
        false_pair_reject = if (scenario == "partial_null") {
          any(fit$table$p_closed[fit$table$contrast != "A - B"] <= alpha)
        } else NA,
        maximum_statistic_difference = max(statistic_difference),
        maximum_p_difference = max(abs(
          fit$intersections$p_perm - direct_intersection_p
        ), abs(fit$table$p_closed - direct_closed)),
        preservation_violations = preservation_violations,
        stringsAsFactors = FALSE
      )
    }
  }
  raw <- do.call(rbind, rows)
  summary <- do.call(rbind, lapply(names(scenarios), function(scenario) {
    values <- raw[raw$scenario == scenario, , drop = FALSE]
    uncertainty <- .binomial_summary(
      sum(values$true_null_reject), nrow(values)
    )
    data.frame(
      scenario = scenario, simulations = nrow(values),
      strong_fwer = uncertainty[["rejection_rate"]],
      monte_carlo_se = uncertainty[["monte_carlo_se"]],
      confidence_low = uncertainty[["confidence_low"]],
      confidence_high = uncertainty[["confidence_high"]],
      false_pair_power = if (scenario == "partial_null") {
        mean(values$false_pair_reject)
      } else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
  rownames(summary) <- NULL
  out <- list(
    summary = summary,
    max_direct_statistic_difference = max(raw$maximum_statistic_difference),
    max_direct_p_difference = max(raw$maximum_p_difference),
    component_preservation_violations = sum(raw$preservation_violations),
    results = if (keep_results) raw else NULL,
    settings = list(
      simulations = simulations, permutations = as.integer(permutations),
      effect = effect, alpha = alpha, seed = seed,
      observations_per_group = per_group,
      permutation_seed_offset = 50000000L
    ),
    proof_boundary = paste(
      "The partial-null audit uses identical complete distributions for A",
      "and B and keeps the non-null C group fixed. It does not validate",
      "centroid equality under unequal distributions."
    ),
    elapsed_seconds = proc.time()[["elapsed"]] - start
  )
  class(out) <- "permanova_strong_pairwise_validation"
  out
}

#' @export
print.permanova_strong_pairwise_validation <- function(x, ...) {
  cat("Strong closed pairwise validation\n")
  print(x$summary, row.names = FALSE)
  cat("Maximum direct statistic difference:",
      x$max_direct_statistic_difference, "\n")
  cat("Component-preservation violations:",
      x$component_preservation_violations, "\n")
  cat("Proof boundary:", x$proof_boundary, "\n")
  invisible(x)
}

#' Benchmark balanced mixed inference
#'
#' @param group_counts Whole-plot counts to benchmark (must be even).
#' @param replicates Observations per whole plot.
#' @param permutations Permutation counts to benchmark.
#' @param seed Reproducibility seed.
#' @return Timing data; results are machine-specific and make no performance claim.
#' @export
benchmark_permanova <- function(group_counts = c(8L, 16L), replicates = 3L,
                                permutations = c(99L, 999L), seed = 1L) {
  group_counts <- as.integer(group_counts)
  permutations <- as.integer(permutations)
  if (any(group_counts < 4L | group_counts %% 2L != 0L)) {
    stop("Every `group_counts` value must be an even integer of at least four.", call. = FALSE)
  }
  replicates <- .positive_integer(replicates, "replicates")
  if (any(permutations < 0L)) stop("`permutations` cannot be negative.", call. = FALSE)
  settings <- expand.grid(groups = group_counts, permutations = permutations)
  rows <- lapply(seq_len(nrow(settings)), function(k) {
    generated <- simulate_balanced_mixed(
      groups_per_treatment = settings$groups[[k]] / 2L,
      replicates = replicates, seed = seed + k
    )
    design <- compile_permanova_design(
      ~ treatment + (1 | site), generated$data, generated$distance
    )
    elapsed <- system.time(suppressWarnings(permanova_mixed(
      design, permutations = settings$permutations[[k]], seed = seed + k
    )))[["elapsed"]]
    data.frame(
      observations = nrow(generated$data), whole_plots = settings$groups[[k]],
      replicates = replicates, permutations = settings$permutations[[k]],
      elapsed_seconds = unname(elapsed), stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
