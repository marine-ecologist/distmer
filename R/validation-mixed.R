# Calibration and exact-agreement validators for the mixed engines added after
# the original support matrix was written. Each engine is pinned against an
# independently computed classical reference, and each reports type-I error and
# power under an exchangeable null and a location shift.
#
# `run_permanova_validation()` calls these with very small simulation counts and
# reads the exact-agreement metric; the rejection-rate summaries are for running
# a validator directly with a realistic number of simulations.

.mixed_validation_summary <- function(raw, simulations, alpha) {
  summary <- do.call(rbind, lapply(split(raw, raw$structure), function(block) {
    do.call(rbind, lapply(split(block, block$scenario), function(part) {
      do.call(rbind, lapply(c("engine", "classical"), function(method) {
        values <- part[[paste0(method, "_p")]]
        values <- values[is.finite(values)]
        rejected <- sum(values <= alpha)
        uncertainty <- .binomial_summary(rejected, max(length(values), 1L))
        data.frame(
          structure = part$structure[[1L]], scenario = part$scenario[[1L]],
          method = method,
          estimand = if (part$scenario[[1L]] == "exchangeable_null") {
            "type_I_error"
          } else "power",
          simulations = length(values), rejected = rejected,
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
  summary
}

.mixed_validation_checks <- function(simulations, permutations, alpha) {
  simulations <- .positive_integer(simulations, "simulations")
  permutations <- .positive_integer(permutations, "permutations")
  if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must lie strictly between zero and one.", call. = FALSE)
  }
  list(simulations = simulations, permutations = permutations)
}

# Shared driver: `step` returns one data frame row per iteration and is handed
# the scenario shift plus a per-iteration seed.
.mixed_validation_run <- function(step, simulations, permutations, effect,
                                  alpha, seed, keep_results, extra = list()) {
  scenarios <- c(exchangeable_null = 0, location_shift = effect)
  rows <- list()
  start <- proc.time()[["elapsed"]]
  for (index in seq_along(scenarios)) {
    scenario <- names(scenarios)[[index]]
    for (iteration in seq_len(simulations)) {
      iteration_seed <- as.integer(seed + index * 100000L + iteration)
      piece <- step(scenarios[[index]], iteration_seed, permutations)
      piece$scenario <- scenario
      rows[[length(rows) + 1L]] <- piece
    }
  }
  raw <- do.call(rbind, rows)
  c(list(
    summary = .mixed_validation_summary(raw, simulations, alpha),
    max_engine_classical_f_difference = max(
      abs(raw$engine_f - raw$classical_f), na.rm = TRUE
    ),
    results = if (keep_results) raw else NULL,
    settings = list(simulations = simulations, permutations = permutations,
                    effect = effect, alpha = alpha, seed = seed),
    elapsed_seconds = proc.time()[["elapsed"]] - start
  ), extra)
}

# ---------------------------------------------------------------------------
# Balanced split-plot / randomized-block / ANCOVA mixed engine
# ---------------------------------------------------------------------------

.mixed_covariate_layout <- function() {
  layout <- expand.grid(rep = 1:2, sub = factor(c("s1", "s2")),
                        unit = factor(seq_len(8L)))
  layout$tr <- factor(ifelse(as.integer(layout$unit) <= 4L, "control", "impact"))
  layout
}

#' Validate the balanced split-plot mixed engine
#'
#' Compares [permanova_mixed_covariate()] with the classical split-plot
#' `aov(y ~ tr * sub + Error(unit))` decomposition: the whole-plot factor is
#' tested in the unit stratum and the sub-plot factor and interaction in the
#' within-unit stratum.
#'
#' @param simulations Iterations per scenario.
#' @param permutations Permutations per fit.
#' @param effect Location shift applied in the alternative scenario.
#' @param unit_sd,residual_sd Simulation standard deviations.
#' @param alpha Nominal level.
#' @param seed Base random seed.
#' @param keep_results Retain the per-iteration rows.
#' @return An object of class `permanova_mixed_covariate_validation`.
#' @export
simulate_mixed_covariate_validation <- function(
    simulations = 200L, permutations = 99L, effect = 1,
    unit_sd = 1, residual_sd = 0.5, alpha = 0.05,
    seed = 20260801L, keep_results = FALSE) {
  checks <- .mixed_validation_checks(simulations, permutations, alpha)
  layout <- .mixed_covariate_layout()
  step <- function(shift, iteration_seed, permutations) {
    set.seed(iteration_seed)
    response <- shift * (layout$tr == "impact") +
      stats::rnorm(nlevels(layout$unit), sd = unit_sd)[layout$unit] +
      stats::rnorm(nrow(layout), sd = residual_sd)
    frame <- layout
    frame$y <- response
    fit <- permanova_mixed_covariate(
      compile_permanova_design(~ tr * sub + (1 | unit), frame, stats::dist(response)),
      permutations = permutations, seed = iteration_seed + 5000L
    )
    reference <- summary(stats::aov(y ~ tr * sub + Error(unit), frame))
    data.frame(
      structure = "balanced split-plot",
      engine_f = fit$table$pseudo_f[fit$table$term == "tr"],
      classical_f = reference[["Error: unit"]][[1L]]["tr", "F value"],
      engine_p = fit$table$p_perm[fit$table$term == "tr"],
      classical_p = reference[["Error: unit"]][[1L]]["tr", "Pr(>F)"],
      stringsAsFactors = FALSE
    )
  }
  out <- .mixed_validation_run(step, checks$simulations, checks$permutations,
                              effect, alpha, seed, keep_results)
  class(out) <- "permanova_mixed_covariate_validation"
  out
}

# ---------------------------------------------------------------------------
# Unbalanced multi-factor split-plot engine
# ---------------------------------------------------------------------------

.unbalanced_splitplot_layout <- function() {
  # Replication is unequal, but every unit size carries the same number of
  # control and impact units. Confounding treatment with unit size would freeze
  # whole size strata under the restricted permutation and make the test
  # conservative for reasons that have nothing to do with the engine.
  sizes <- c(2L, 3L, 4L, 5L)
  rows <- list()
  unit <- 0L
  for (size in sizes) {
    for (treatment in c("control", "impact")) {
      for (replicate in 1:2) {
        unit <- unit + 1L
        for (sub in c("s1", "s2")) {
          rows[[length(rows) + 1L]] <- data.frame(
            unit = sprintf("u%02d", unit), sub = sub,
            rep = seq_len(size), tr = treatment, stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  layout <- do.call(rbind, rows)
  layout$unit <- factor(layout$unit)
  layout$sub <- factor(layout$sub)
  layout$tr <- factor(layout$tr)
  layout
}

#' Validate the unbalanced split-plot mixed engine
#'
#' Compares [permanova_unbalanced_splitplot()] with an `lm` on equally weighted
#' unit centroids for the whole-plot factor, which is the estimand the engine
#' declares under unequal replication.
#'
#' @inheritParams simulate_mixed_covariate_validation
#' @return An object of class `permanova_unbalanced_splitplot_validation`.
#' @export
simulate_unbalanced_splitplot_validation <- function(
    simulations = 200L, permutations = 99L, effect = 1,
    unit_sd = 1, residual_sd = 0.5, alpha = 0.05,
    seed = 20260802L, keep_results = FALSE) {
  checks <- .mixed_validation_checks(simulations, permutations, alpha)
  layout <- .unbalanced_splitplot_layout()
  step <- function(shift, iteration_seed, permutations) {
    set.seed(iteration_seed)
    response <- shift * (layout$tr == "impact") +
      stats::rnorm(nlevels(layout$unit), sd = unit_sd)[layout$unit] +
      stats::rnorm(nrow(layout), sd = residual_sd)
    frame <- layout
    frame$y <- response
    fit <- permanova_unbalanced_splitplot(
      compile_permanova_design(~ tr + sub + (1 | unit), frame, stats::dist(response)),
      permutations = permutations, seed = iteration_seed + 5000L
    )
    centroids <- stats::aggregate(y ~ unit + tr, frame, mean)
    reference <- stats::anova(stats::lm(y ~ tr, centroids))
    data.frame(
      structure = "unbalanced split-plot",
      engine_f = fit$table$pseudo_f[fit$table$term == "tr"],
      classical_f = reference["tr", "F value"],
      engine_p = fit$table$p_perm[fit$table$term == "tr"],
      classical_p = reference["tr", "Pr(>F)"],
      stringsAsFactors = FALSE
    )
  }
  out <- .mixed_validation_run(step, checks$simulations, checks$permutations,
                              effect, alpha, seed, keep_results)
  class(out) <- "permanova_unbalanced_splitplot_validation"
  out
}

# ---------------------------------------------------------------------------
# Covariate-adjusted balanced nested engine
# ---------------------------------------------------------------------------

.nested_covariate_layout <- function() {
  layout <- expand.grid(rep = 1:3, site = factor(1:2), reg = factor(1:6))
  layout$tr <- factor(ifelse(as.integer(layout$reg) <= 3L, "control", "impact"))
  layout
}

#' Validate the covariate-adjusted balanced nested engine
#'
#' [permanova_nested_covariate()] charges the covariate to the observation
#' residual and tests each factor for the sum of squares it adds *over* the
#' covariate, against the exact nested stratum of the factor-only model. This
#' validator reconstructs both pieces independently: the numerator as
#' `SS(factor | covariate)` from an observation-level `lm`, and the denominator
#' as the corresponding stratum mean square of `aov(y ~ tr + Error(reg/site))`.
#'
#' @inheritParams simulate_mixed_covariate_validation
#' @param region_sd,site_sd Simulation standard deviations for the hierarchy.
#' @param residual_sd Observation-level simulation standard deviation.
#' @param slope Covariate slope used to generate the response.
#' @return An object of class `permanova_nested_covariate_validation`.
#' @export
simulate_nested_covariate_validation <- function(
    simulations = 200L, permutations = 99L, effect = 1,
    region_sd = 1, site_sd = 0.5, residual_sd = 0.4, slope = 0.7,
    alpha = 0.05, seed = 20260803L, keep_results = FALSE) {
  checks <- .mixed_validation_checks(simulations, permutations, alpha)
  layout <- .nested_covariate_layout()
  site_label <- interaction(layout$reg, layout$site, drop = TRUE)
  step <- function(shift, iteration_seed, permutations) {
    set.seed(iteration_seed)
    frame <- layout
    frame$x <- stats::rnorm(nrow(layout))
    frame$y <- shift * (layout$tr == "impact") + slope * frame$x +
      stats::rnorm(nlevels(layout$reg), sd = region_sd)[layout$reg] +
      stats::rnorm(nlevels(site_label), sd = site_sd)[site_label] +
      stats::rnorm(nrow(layout), sd = residual_sd)
    fit <- permanova_nested_covariate(
      compile_permanova_design(~ x + tr + (1 | reg / site), frame,
                               stats::dist(frame$y)),
      permutations = permutations, seed = iteration_seed + 5000L
    )
    # numerator: sum of squares the factor adds over the covariate
    added <- stats::deviance(stats::lm(y ~ x, frame)) -
      stats::deviance(stats::lm(y ~ x + tr, frame))
    # denominator: region stratum of the factor-only nested model
    strata <- summary(stats::aov(y ~ tr + Error(reg / site), frame))
    region <- strata[["Error: reg"]][[1L]]
    denominator_ms <- region["Residuals", "Mean Sq"]
    denominator_df <- region["Residuals", "Df"]
    classical_f <- added / denominator_ms
    data.frame(
      structure = "balanced nested with covariate",
      engine_f = fit$table$pseudo_f[fit$table$term == "tr"],
      classical_f = classical_f,
      engine_p = fit$table$p_perm[fit$table$term == "tr"],
      classical_p = stats::pf(classical_f, 1L, denominator_df,
                              lower.tail = FALSE),
      stringsAsFactors = FALSE
    )
  }
  out <- .mixed_validation_run(step, checks$simulations, checks$permutations,
                              effect, alpha, seed, keep_results)
  class(out) <- "permanova_nested_covariate_validation"
  out
}

# ---------------------------------------------------------------------------
# Asymmetric nested mixed engine (three error strata)
# ---------------------------------------------------------------------------

.asymmetric_mixed_layout <- function() {
  layout <- expand.grid(rep = 1:3, DE = factor(c("shallow", "deep")),
                        ISLAND = factor(paste0("i", seq_len(8L))))
  layout$EC <- factor(ifelse(as.integer(layout$ISLAND) <= 4L, "e1", "e2"))
  layout
}

#' Validate the asymmetric nested mixed engine
#'
#' Compares [permanova_asymmetric_mixed()] with
#' `aov(y ~ EC * DE + Error(ISLAND / DE))`: the whole-plot factor is tested in
#' the unit stratum, and the sub-plot factor and its interaction in the
#' unit-by-sub-plot stratum.
#'
#' @inheritParams simulate_mixed_covariate_validation
#' @param island_sd,island_depth_sd Simulation standard deviations for the two
#'   random strata.
#' @param residual_sd Observation-level simulation standard deviation.
#' @return An object of class `permanova_asymmetric_mixed_validation`.
#' @export
simulate_asymmetric_mixed_validation <- function(
    simulations = 200L, permutations = 99L, effect = 1,
    island_sd = 1, island_depth_sd = 0.6, residual_sd = 0.4,
    alpha = 0.05, seed = 20260804L, keep_results = FALSE) {
  checks <- .mixed_validation_checks(simulations, permutations, alpha)
  layout <- .asymmetric_mixed_layout()
  cell <- interaction(layout$ISLAND, layout$DE, drop = TRUE)
  step <- function(shift, iteration_seed, permutations) {
    set.seed(iteration_seed)
    frame <- layout
    frame$y <- shift * (layout$EC == "e2") +
      stats::rnorm(nlevels(layout$ISLAND), sd = island_sd)[layout$ISLAND] +
      stats::rnorm(nlevels(cell), sd = island_depth_sd)[cell] +
      stats::rnorm(nrow(layout), sd = residual_sd)
    fit <- permanova_asymmetric_mixed(
      compile_permanova_design(~ EC * DE + (1 | ISLAND) + (1 | ISLAND:DE),
                               frame, stats::dist(frame$y)),
      permutations = permutations, seed = iteration_seed + 5000L
    )
    reference <- summary(stats::aov(y ~ EC * DE + Error(ISLAND / DE), frame))
    island <- reference[["Error: ISLAND"]][[1L]]
    data.frame(
      structure = "asymmetric nested mixed",
      engine_f = fit$table$pseudo_f[fit$table$term == "EC"],
      classical_f = island["EC", "F value"],
      engine_p = fit$table$p_perm[fit$table$term == "EC"],
      classical_p = island["EC", "Pr(>F)"],
      stringsAsFactors = FALSE
    )
  }
  out <- .mixed_validation_run(step, checks$simulations, checks$permutations,
                              effect, alpha, seed, keep_results)
  class(out) <- "permanova_asymmetric_mixed_validation"
  out
}

# ---------------------------------------------------------------------------
# Continuous random-slope engine
# ---------------------------------------------------------------------------

#' Validate the continuous random-slope engine
#'
#' [permanova_random_slope_continuous()] tests the equally weighted mean
#' within-group slope of a numeric covariate against the variation among group
#' slopes. For a Euclidean response that is exactly the one-sample test on the
#' per-group ordinary least squares slopes, so the reference here is
#' `t.test(slopes)` squared.
#'
#' @inheritParams simulate_mixed_covariate_validation
#' @param groups Number of random groups.
#' @param occasions Covariate values observed within every group.
#' @param slope_sd Standard deviation of the group slopes.
#' @param residual_sd Observation-level simulation standard deviation.
#' @return An object of class `permanova_random_slope_continuous_validation`.
#' @export
simulate_random_slope_continuous_validation <- function(
    simulations = 200L, permutations = 99L, effect = 0.8,
    groups = 8L, occasions = 5L, slope_sd = 0.5, residual_sd = 0.3,
    alpha = 0.05, seed = 20260805L, keep_results = FALSE) {
  checks <- .mixed_validation_checks(simulations, permutations, alpha)
  groups <- .positive_integer(groups, "groups")
  occasions <- .positive_integer(occasions, "occasions")
  group <- factor(rep(seq_len(groups), each = occasions))
  covariate <- rep(seq_len(occasions), times = groups)
  step <- function(shift, iteration_seed, permutations) {
    set.seed(iteration_seed)
    slopes <- stats::rnorm(groups, mean = shift, sd = slope_sd)
    frame <- data.frame(g = group, x = covariate)
    frame$y <- slopes[group] * covariate +
      stats::rnorm(length(covariate), sd = residual_sd)
    fit <- permanova_random_slope_continuous(
      compile_permanova_design(~ x + (1 + x | g), frame, stats::dist(frame$y)),
      permutations = permutations, seed = iteration_seed + 5000L
    )
    observed_slopes <- vapply(split(seq_len(nrow(frame)), frame$g), function(i) {
      unname(stats::coef(stats::lm(frame$y[i] ~ frame$x[i]))[[2L]])
    }, numeric(1))
    reference <- stats::t.test(observed_slopes)
    data.frame(
      structure = "continuous random slope",
      engine_f = fit$test$statistic,
      classical_f = unname(reference$statistic^2),
      engine_p = fit$test$p_value,
      classical_p = reference$p.value,
      stringsAsFactors = FALSE
    )
  }
  out <- .mixed_validation_run(step, checks$simulations, checks$permutations,
                              effect, alpha, seed, keep_results)
  class(out) <- "permanova_random_slope_continuous_validation"
  out
}
