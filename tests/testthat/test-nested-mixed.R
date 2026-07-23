.nested_example <- function() {
  dat <- expand.grid(
    replicate = factor(1:3), site = factor(1:2), region = factor(1:6)
  )
  dat$treatment <- factor(
    ifelse(as.integer(dat$region) <= 3, "control", "impact")
  )
  region_effect <- c(-1.2, 0.1, 1.0, -0.8, 0.3, 1.3)
  site_effect <- c(-0.4, 0.4, 0.2, -0.2, -0.1, 0.1,
                   0.3, -0.3, -0.2, 0.2, 0.5, -0.5)
  nested_site <- interaction(dat$region, dat$site, drop = TRUE)
  dat$y <- 1.5 * (dat$treatment == "impact") +
    region_effect[dat$region] + site_effect[nested_site] +
    rep(c(-0.1, 0, 0.1), 12)
  dat
}

.nested_multilevel_example <- function(regions = 6L, replicates = 3L) {
  dat <- expand.grid(
    replicate = factor(seq_len(replicates)), site = factor(1:2),
    region = factor(seq_len(regions))
  )
  dat$region_treatment <- factor(
    ifelse(as.integer(dat$region) <= regions / 2, "control", "impact")
  )
  dat$site_treatment <- dat$site
  dat$observation_treatment <- dat$replicate
  nested_site <- interaction(dat$region, dat$site, drop = TRUE)
  set.seed(81 + regions + replicates)
  observation_effect <- seq(-0.4, 0.4, length.out = replicates)
  dat$y <- 1.2 * (dat$region_treatment == "impact") +
    0.7 * (dat$site_treatment == "2") +
    observation_effect[dat$observation_treatment] +
    rnorm(regions)[dat$region] +
    rnorm(2L * regions, sd = 0.5)[nested_site] +
    rnorm(nrow(dat), sd = 0.2)
  dat
}

test_that("nested EMS matches every random component exactly", {
  dat <- .nested_example()
  design <- compile_permanova_design(
    ~ treatment + (1 | region/site), dat, dist(dat$y)
  )
  ems <- derive_nested_ems(design)$treatment

  expect_equal(design$random_term_table$origin, c("main", "nesting"))
  expect_equal(ems$numerator_under_null,
               c(residual = 1, region = 6, `region:site` = 3), tolerance = 1e-10)
  expect_equal(ems$denominator, ems$numerator_under_null, tolerance = 1e-10)
  expect_equal(ems$denominator_term, "region within fixed model")
  expect_true(ems$exact_denominator)
})

test_that("nested Euclidean pseudo-F equals ANOVA on region means", {
  dat <- .nested_example()
  design <- compile_permanova_design(
    ~ treatment + (1 | region/site), dat, dist(dat$y)
  )
  fit <- suppressWarnings(permanova_nested(design, permutations = 0, truncate = TRUE))
  region_means <- aggregate(y ~ region + treatment, dat, mean)
  reference <- summary(aov(y ~ treatment, region_means))[[1]][["F value"]][1]

  expect_s3_class(fit, "permanova_nested")
  expect_equal(fit$tests$treatment$statistic, unname(reference), tolerance = 1e-10)
  expect_equal(fit$tests$treatment$exchangeable_units, "region")
  expect_equal(fit$hierarchy, c("region", "region:site"))
  expect_named(fit$variance_components$unrestricted,
               c("residual", "region", "region:site"))
  expect_equal(
    as.vector(fit$variance_components$coefficient_matrix %*%
      fit$variance_components$unrestricted),
    unname(fit$variance_components$observed_mean_squares), tolerance = 1e-10
  )
  expect_true(all(fit$variance_components$truncated >= 0))
})

test_that("fixed terms are assigned to region, site, and observation strata", {
  dat <- .nested_multilevel_example()
  design <- compile_permanova_design(
    ~ region_treatment + site_treatment + observation_treatment +
      (1 | region/site),
    dat, dist(dat$y)
  )
  fit <- permanova_nested(design, permutations = 0)

  expect_equal(
    fit$assignment_levels,
    c(
      region_treatment = "region", site_treatment = "region:site",
      observation_treatment = "observation"
    )
  )
  expect_equal(
    fit$ems$region_treatment$numerator_under_null,
    c(residual = 1, region = 6, `region:site` = 3), tolerance = 1e-10
  )
  expect_equal(
    fit$ems$site_treatment$numerator_under_null,
    c(residual = 1, region = 0, `region:site` = 3), tolerance = 1e-10
  )
  expect_equal(
    fit$ems$observation_treatment$numerator_under_null,
    c(residual = 1, region = 0, `region:site` = 0), tolerance = 1e-10
  )
  expect_equal(
    vapply(fit$tests, `[[`, character(1), "parent_level"),
    c(
      region_treatment = "grand mean", site_treatment = "region",
      observation_treatment = "region:site"
    )
  )
})

test_that("level-specific pseudo-F statistics equal classical nested ANOVA", {
  dat <- .nested_multilevel_example()
  design <- compile_permanova_design(
    ~ region_treatment + site_treatment + observation_treatment +
      (1 | region/site),
    dat, dist(dat$y)
  )
  fit <- permanova_nested(design, permutations = 0)

  region_means <- aggregate(y ~ region + region_treatment, dat, mean)
  site_means <- aggregate(y ~ region + site + site_treatment, dat, mean)
  nested_site <- interaction(dat$region, dat$site, drop = TRUE)
  region_reference <- summary(aov(
    y ~ region_treatment, region_means
  ))[[1]]
  site_reference <- summary(aov(
    y ~ region + site_treatment, site_means
  ))[[1]]
  observation_reference <- summary(aov(
    dat$y ~ nested_site + dat$observation_treatment
  ))[[1]]

  expect_equal(
    fit$tests$region_treatment$statistic,
    unname(region_reference["region_treatment", "F value"]),
    tolerance = 1e-10
  )
  expect_equal(
    fit$tests$site_treatment$statistic,
    unname(site_reference["site_treatment", "F value"]),
    tolerance = 1e-10
  )
  expect_equal(
    fit$tests$observation_treatment$statistic,
    unname(observation_reference["dat$observation_treatment", "F value"]),
    tolerance = 1e-10
  )
})

test_that("balanced cross-level fixed interactions inherit the lower stratum", {
  dat <- .nested_multilevel_example()
  contrasts(dat$region_treatment) <- contr.sum(2)
  contrasts(dat$site_treatment) <- contr.sum(2)
  design <- compile_permanova_design(
    ~ region_treatment * site_treatment + (1 | region/site),
    dat, dist(dat$y)
  )
  fit <- permanova_nested(design, permutations = 0)

  expect_equal(
    fit$assignment_levels,
    c(
      region_treatment = "region", site_treatment = "region:site",
      `region_treatment:site_treatment` = "region:site"
    )
  )
  expect_equal(
    fit$ems$`region_treatment:site_treatment`$numerator_under_null,
    c(residual = 1, region = 0, `region:site` = 3),
    tolerance = 1e-10
  )
})

test_that("cross-level contrast rejection names the explicit remedy", {
  dat <- .nested_multilevel_example()
  contrasts(dat$region_treatment) <- contr.treatment(2)
  contrasts(dat$site_treatment) <- contr.treatment(2)
  design <- compile_permanova_design(
    ~ region_treatment * site_treatment + (1 | region/site),
    dat, dist(dat$y)
  )

  expect_error(
    permanova_nested(design, permutations = 0),
    "explicit orthogonal contrasts such as `contr.sum`"
  )
})

test_that("small nested hierarchy enumerates top units exactly", {
  dat <- .nested_example()
  design <- compile_permanova_design(
    ~ treatment + (1 | region/site), dat, dist(dat$y)
  )
  fit <- suppressWarnings(permanova_nested(
    design, permutations = 719, seed = 3, keep_permutations = TRUE
  ))

  expect_true(fit$permutation$exact)
  expect_equal(fit$permutation$attainable, factorial(6))
  expect_equal(fit$permutation$used, factorial(6))
  expect_equal(fit$tests$treatment$p_value * factorial(6) %% 1,
               0, tolerance = 1e-12)
  expect_equal(length(fit$tests$treatment$permutation_statistics), factorial(6))
})

test_that("site and observation permutations enumerate within parents", {
  dat <- .nested_multilevel_example(regions = 3, replicates = 2)
  site_design <- compile_permanova_design(
    ~ site_treatment + (1 | region/site), dat, dist(dat$y)
  )
  site_fit <- suppressWarnings(permanova_nested(
    site_design, permutations = "all", keep_permutations = TRUE
  ))
  expect_equal(site_fit$tests$site_treatment$assignment_level, "region:site")
  expect_equal(site_fit$tests$site_treatment$parent_level, "region")
  expect_equal(site_fit$tests$site_treatment$attainable_permutations, 2 ^ 3)
  expect_equal(site_fit$tests$site_treatment$permutations_used, 2 ^ 3)
  expect_equal(site_fit$tests$site_treatment$unique_permutations_used, 2 ^ 3)
  expect_equal(site_fit$tests$site_treatment$p_value * 2 ^ 3 %% 1,
               0, tolerance = 1e-12)

  observation_design <- compile_permanova_design(
    ~ observation_treatment + (1 | region/site), dat, dist(dat$y)
  )
  observation_fit <- suppressWarnings(permanova_nested(
    observation_design, permutations = "all", keep_permutations = TRUE
  ))
  expect_equal(
    observation_fit$tests$observation_treatment$assignment_level,
    "observation"
  )
  expect_equal(
    observation_fit$tests$observation_treatment$parent_level,
    "region:site"
  )
  expect_equal(
    observation_fit$tests$observation_treatment$attainable_permutations,
    2 ^ 6
  )
  expect_equal(
    observation_fit$tests$observation_treatment$unique_permutations_used,
    2 ^ 6
  )

  site_structure <- distmer:::.validate_nested_design(site_design)
  site_plan <- distmer:::.nested_permutation_plan(
    site_structure, level = 2L, permutations = "all", max_exact = 1000L
  )
  for (units in site_plan$units_by_parent) {
    expect_true(all(apply(site_plan$orders[, units, drop = FALSE], 1, function(x) {
      identical(sort(as.integer(x)), sort(as.integer(units)))
    })))
  }
})

test_that("unsupported nested hypotheses and imbalance fail explicitly", {
  dat <- .nested_example()
  dat$ambiguous <- as.numeric(dat$region) + as.numeric(dat$replicate)
  ambiguous <- compile_permanova_design(
    ~ ambiguous + (1 | region/site), dat, dist(dat$y)
  )
  expect_error(permanova_nested(ambiguous, permutations = 0),
               "spans multiple hierarchy levels")

  unbalanced <- dat[-1L, ]
  unbalanced_design <- compile_permanova_design(
    ~ treatment + (1 | region/site), unbalanced, dist(unbalanced$y)
  )
  expect_error(permanova_nested(unbalanced_design, permutations = 0),
               "equal observations")

  crossed <- compile_permanova_design(
    ~ treatment + (1 | region) + (1 | site), dat, dist(dat$y)
  )
  expect_error(permanova_nested(crossed, permutations = 0), "/`-generated hierarchy")

  saturated <- .nested_multilevel_example(regions = 2, replicates = 2)
  saturated_design <- compile_permanova_design(
    ~ region_treatment + (1 | region/site), saturated, dist(saturated$y)
  )
  expect_error(permanova_nested(saturated_design, permutations = 0),
               "No residual degrees of freedom")

  branched <- compile_permanova_design(
    ~ treatment + (1 | region/site) + (1 | region/replicate),
    dat, dist(dat$y)
  )
  expect_error(permanova_nested(branched, permutations = 0),
               "contiguous containment chain")
})

test_that("three-level nested hierarchies compile and estimate components", {
  dat <- expand.grid(
    replicate = factor(1:2), quadrat = factor(1:2),
    site = factor(1:2), region = factor(1:6)
  )
  dat$treatment <- factor(ifelse(as.integer(dat$region) <= 3, "a", "b"))
  set.seed(31)
  dat$y <- rnorm(nrow(dat))
  design <- compile_permanova_design(
    ~ treatment + (1 | region/site/quadrat), dat, dist(dat$y)
  )
  fit <- permanova_nested(design, permutations = 0)

  expect_equal(fit$hierarchy,
               c("region", "region:site", "region:site:quadrat"))
  expect_named(fit$variance_components$unrestricted,
               c("residual", "region", "region:site", "region:site:quadrat"))
})
