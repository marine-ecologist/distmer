.hybrid_example <- function(regions = 4L, sites = 2L, replicates = 2L,
                            treatments = 2L) {
  dat <- expand.grid(
    replicate = factor(seq_len(replicates)),
    treatment = factor(seq_len(treatments)),
    site = factor(seq_len(sites)), region = factor(seq_len(regions))
  )
  region_site <- interaction(dat$region, dat$site, drop = TRUE)
  treatment_region <- interaction(dat$treatment, dat$region, drop = TRUE)
  treatment_site <- interaction(
    dat$treatment, dat$region, dat$site, drop = TRUE
  )
  set.seed(501 + regions + sites + treatments)
  dat$y <-
    rnorm(regions, sd = 1)[dat$region] +
    rnorm(regions * sites, sd = 0.6)[region_site] +
    rnorm(treatments * regions, sd = 0.8)[treatment_region] +
    rnorm(treatments * regions * sites, sd = 0.4)[treatment_site] +
    rnorm(nrow(dat), sd = 0.2)
  dat
}

.hybrid_formula <- function() {
  ~ treatment + (1 | region/site) +
    (1 | treatment:region) + (1 | treatment:region:site)
}

.hybrid_profile_f <- function(response, data) {
  means <- aggregate(
    response,
    list(treatment = data$treatment, region = data$region), mean
  )
  table <- summary(aov(x ~ treatment * region, means))[[1]]
  rows <- trimws(rownames(table))
  unname(table[match("treatment", rows), "Mean Sq"] /
           table[match("treatment:region", rows), "Mean Sq"])
}

test_that("hybrid EMS gives the exact treatment-by-region denominator", {
  dat <- .hybrid_example()
  design <- compile_permanova_design(
    .hybrid_formula(), dat, dist(dat$y)
  )
  ems <- derive_hybrid_ems(design)

  expect_equal(
    ems$numerator_under_null,
    c(
      residual = 1, region = 0, `region:site` = 0,
      `treatment:region` = 4, `treatment:region:site` = 2
    ),
    tolerance = 1e-10
  )
  expect_equal(
    ems$numerator_under_null, ems$denominator, tolerance = 1e-10
  )
  expect_equal(ems$denominator_term, "treatment:region")
  expect_true(ems$exact_denominator)
})

test_that("hybrid Euclidean F equals ANOVA on treatment-region profiles", {
  dat <- .hybrid_example(regions = 6L, sites = 3L, replicates = 2L)
  design <- compile_permanova_design(
    .hybrid_formula(), dat, dist(dat$y)
  )
  fit <- permanova_hybrid(design, permutations = 0, truncate = TRUE)
  reference <- .hybrid_profile_f(dat$y, dat)

  expect_s3_class(fit, "permanova_hybrid")
  expect_equal(fit$test$statistic, reference, tolerance = 1e-10)
  expect_equal(fit$test$denominator$expression, "MS[treatment:region]")
  expect_equal(fit$test$exchangeable_units, "treatment:region")
  expect_equal(fit$test$permutation_strata, "region")
  expect_true(all(fit$variance_components$truncated >= 0))
})

test_that("hybrid exact permutations match every direct profile ANOVA", {
  dat <- .hybrid_example(regions = 3L)
  design <- compile_permanova_design(
    .hybrid_formula(), dat, dist(dat$y)
  )
  fit <- suppressWarnings(permanova_hybrid(
    design, permutations = "all", keep_permutations = TRUE
  ))
  structure <- distmer:::.validate_hybrid_design(design)
  plan <- distmer:::.hybrid_permutation_plan(
    structure, "all", 100L
  )
  direct <- vapply(seq_len(nrow(plan$orders)), function(row) {
    index <- distmer:::.block_permutation_index(
      structure$profile_blocks, plan$orders[row, ]
    )
    .hybrid_profile_f(dat$y[index], dat)
  }, numeric(1))

  expect_true(fit$test$exact_permutation)
  expect_equal(fit$test$attainable_permutations, 2 ^ 3)
  expect_equal(fit$test$unique_permutations_used, 2 ^ 3)
  expect_equal(fit$test$permutation_statistics, direct, tolerance = 1e-10)
  expect_equal(fit$test$p_value * 2 ^ 3 %% 1, 0, tolerance = 1e-12)
})

test_that("hybrid permutations keep profiles and descendants intact", {
  dat <- .hybrid_example(regions = 3L, sites = 3L, replicates = 2L)
  design <- compile_permanova_design(
    .hybrid_formula(), dat, dist(dat$y)
  )
  structure <- distmer:::.validate_hybrid_design(design)
  plan <- distmer:::.hybrid_permutation_plan(
    structure, "all", 100L
  )

  for (row in seq_len(nrow(plan$orders))) {
    for (region_profiles in structure$profiles_by_region) {
      expect_setequal(plan$orders[row, region_profiles], region_profiles)
    }
    index <- distmer:::.block_permutation_index(
      structure$profile_blocks, plan$orders[row, ]
    )
    for (target in seq_along(structure$profile_blocks)) {
      source <- plan$orders[row, target]
      expect_equal(
        index[structure$profile_blocks[[target]]],
        structure$profile_blocks[[source]]
      )
    }
  }
})

test_that("hybrid variance components reconstruct every selected mean square", {
  dat <- .hybrid_example(regions = 5L, sites = 3L)
  design <- compile_permanova_design(
    .hybrid_formula(), dat, dist(dat$y)
  )
  components <- estimate_hybrid_variance_components(design)

  expect_named(
    components$unrestricted,
    c(
      "residual", "region", "region:site", "treatment:region",
      "treatment:region:site"
    )
  )
  expect_equal(
    as.vector(components$coefficient_matrix %*% components$unrestricted),
    unname(components$observed_mean_squares), tolerance = 1e-10
  )
})

test_that("hybrid Monte Carlo permutations are seed reproducible", {
  dat <- .hybrid_example(
    regions = 4L, sites = 2L, replicates = 2L, treatments = 3L
  )
  design <- compile_permanova_design(
    .hybrid_formula(), dat, dist(dat$y)
  )
  first <- permanova_hybrid(
    design, permutations = 19, seed = 91, keep_permutations = TRUE
  )
  second <- permanova_hybrid(
    design, permutations = 19, seed = 91, keep_permutations = TRUE
  )

  expect_false(first$test$exact_permutation)
  expect_equal(first$test$attainable_permutations, factorial(3) ^ 4)
  expect_equal(
    first$test$permutation_statistics,
    second$test$permutation_statistics
  )
  expect_equal(first$test$p_value, second$test$p_value)
})

test_that("hybrid non-Euclidean distances retain empirical inference", {
  dat <- .hybrid_example(regions = 3L)
  distance <- as.matrix(dist(cbind(seq_len(nrow(dat)), dat$y)))
  distance[1L, 2L] <- distance[2L, 1L] <- 4 * max(distance)
  design <- compile_permanova_design(.hybrid_formula(), dat, distance)
  fit <- suppressWarnings(permanova_hybrid(
    design, permutations = "all"
  ))

  expect_false(design$spectral$euclidean)
  expect_true(is.finite(fit$test$p_value))
  expect_match(
    paste(fit$test$warnings, collapse = " "),
    "complete non-Euclidean"
  )
  expect_equal(
    fit$design$spectral$eigenvalues,
    design$spectral$eigenvalues
  )
})

test_that("unsupported hybrid structures fail before inference", {
  dat <- .hybrid_example()
  distance <- dist(dat$y)

  missing_region_interaction <- compile_permanova_design(
    ~ treatment + (1 | region/site) + (1 | treatment:region:site),
    dat, distance
  )
  expect_error(
    permanova_hybrid(missing_region_interaction, permutations = 0),
    "requires region, region/site"
  )

  missing_site_interaction <- compile_permanova_design(
    ~ treatment + (1 | region/site) + (1 | treatment:region),
    dat, distance
  )
  expect_error(
    permanova_hybrid(missing_site_interaction, permutations = 0),
    "requires region, region/site"
  )

  one_site <- subset(dat, site == "1")
  one_site <- droplevels(one_site)
  one_site_design <- compile_permanova_design(
    .hybrid_formula(), one_site, dist(one_site$y)
  )
  expect_error(
    permanova_hybrid(one_site_design, permutations = 0),
    "at least two sites"
  )

  unbalanced <- dat[-1L, ]
  unbalanced_design <- compile_permanova_design(
    .hybrid_formula(), unbalanced, dist(unbalanced$y)
  )
  expect_error(
    permanova_hybrid(unbalanced_design, permutations = 0),
    "equal replication"
  )

  continuous <- dat
  continuous$treatment <- as.numeric(continuous$treatment)
  continuous_design <- compile_permanova_design(
    .hybrid_formula(), continuous, distance
  )
  expect_error(
    permanova_hybrid(continuous_design, permutations = 0),
    "categorical fixed treatment"
  )
})
