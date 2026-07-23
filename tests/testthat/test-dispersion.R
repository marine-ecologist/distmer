.dispersion_mixed_example <- function() {
  dat <- expand.grid(replicate = factor(1:3), cluster = factor(1:8))
  dat$treatment <- factor(
    ifelse(as.integer(dat$cluster) <= 4, "control", "impact")
  )
  set.seed(401)
  cluster_effect <- matrix(rnorm(8 * 3, sd = 0.4), 8, 3)
  scale <- ifelse(dat$treatment == "impact", 1.8, 0.7)
  response <- cluster_effect[dat$cluster, ] +
    matrix(rnorm(nrow(dat) * 3), ncol = 3) * scale
  list(data = dat, response = response)
}

.dispersion_nested_example <- function() {
  dat <- expand.grid(
    replicate = factor(1:2), site = factor(1:2), region = factor(1:3)
  )
  dat$treatment <- dat$site
  set.seed(402)
  response <- matrix(rnorm(nrow(dat) * 3), ncol = 3) *
    ifelse(dat$treatment == "1", 0.7, 1.5)
  list(data = dat, response = response)
}

.dispersion_crossed_example <- function() {
  dat <- expand.grid(
    replicate = factor(1:2), treatment = factor(c("a", "b")),
    block = factor(1:3)
  )
  set.seed(403)
  response <- matrix(rnorm(nrow(dat) * 4), ncol = 4) *
    ifelse(dat$treatment == "a", 0.8, 1.4)
  list(data = dat, response = response)
}

test_that("one-way dispersion geometry and permutation match vegan", {
  skip_if_not_installed("vegan")
  data(dune, dune.env, package = "vegan")
  complete <- dune.env$Management %in% c("NM", "SF")
  data <- droplevels(dune.env[complete, , drop = FALSE])
  distance <- as.dist(as.matrix(vegan::vegdist(dune))[complete, complete])
  design <- compile_permanova_design(~ Management, data, distance)

  for (center in c("centroid", "median")) {
    reference <- vegan::betadisper(
      distance, data$Management, type = center
    )
    fit <- permanova_dispersion(
      distance ~ Management, data, type = center, permutations = 19,
      seed = 71, keep_permutations = TRUE
    )
    expect_equal(fit$distances, reference$distances, tolerance = 1e-9)
    expect_equal(
      fit$test$statistic,
      unname(anova(reference)["Groups", "F value"]),
      tolerance = 1e-9
    )

    geometry <- distmer:::.dispersion_geometry(
      distance, data$Management, type = center
    )
    response_design <- distmer:::.dispersion_response_design(
      design, geometry$distances, 1e-10
    )
    set.seed(71)
    setup <- distmer:::.dispersion_plan(
      design, response_design, 19L, 100000L, 1e-10
    )
    vegan_permutation <- vegan::permutest(
      reference, permutations = setup$plan$orders
    )
    expect_equal(
      fit$test$permutation_statistics,
      unname(vegan_permutation$perm),
      tolerance = 1e-9
    )
    expect_equal(
      fit$test$p_value,
      unname(vegan_permutation$tab["Groups", "Pr(>F)"]),
      tolerance = 1e-12
    )
  }
})

test_that("dispersion accepts formula, one-sided, and compiled-design inputs", {
  skip_if_not_installed("vegan")
  data(dune, dune.env, package = "vegan")
  selected <- dune.env$Management %in% c("NM", "SF")
  env <- droplevels(dune.env[selected, , drop = FALSE])
  bray <- vegan::vegdist(dune[selected, ], method = "bray")
  design <- compile_permanova_design(~ Management, env, bray)

  direct <- permanova_dispersion(
    bray ~ Management, env, permutations = 0
  )
  one_sided <- permanova_dispersion(
    ~ Management, env, distance = bray, permutations = 0
  )
  legacy <- permanova_dispersion(
    design, "Management", permutations = 0
  )
  named_legacy <- permanova_dispersion(
    design = design, term = "Management", permutations = 0
  )

  expect_s3_class(direct, "permanova_dispersion")
  expect_equal(direct$term, "Management")
  expect_equal(direct$test$statistic, one_sided$test$statistic)
  expect_equal(direct$test$statistic, legacy$test$statistic)
  expect_equal(direct$test$statistic, named_legacy$test$statistic)
  expect_equal(direct$design$distance, as.matrix(bray))
  expect_error(
    permanova_dispersion(
      bray ~ Management, env, distance = bray, permutations = 0
    ),
    "Supply the dissimilarity once"
  )
})

test_that("small one-way dispersion space is exhaustive", {
  dat <- data.frame(
    treatment = factor(rep(c("a", "b"), each = 3))
  )
  response <- cbind(
    c(-1.2, -0.1, 0.8, 2.0, 2.4, 4.1),
    c(0.1, 0.5, -0.4, 0.2, 1.3, -0.8)
  )
  design <- compile_permanova_design(
    ~ treatment, dat, dist(response)
  )
  fit <- suppressWarnings(permanova_dispersion(
    design, "treatment", type = "centroid", permutations = "all",
    keep_permutations = TRUE, max_exact = 1000L
  ))

  residuals <- residuals(lm(fit$distances ~ dat$treatment))
  orders <- distmer:::.all_orders(nrow(dat))
  direct <- apply(orders, 1L, function(index) {
    unname(summary(aov(residuals[index] ~ dat$treatment))[[1]][1, "F value"])
  })
  expect_true(fit$test$exact_permutation)
  expect_equal(fit$test$attainable_permutations, factorial(6))
  expect_equal(fit$test$permutation_statistics, direct, tolerance = 1e-10)
  expect_equal(fit$test$p_value * factorial(6) %% 1, 0, tolerance = 1e-12)
})

test_that("random-intercept dispersion uses the whole-plot denominator", {
  example <- .dispersion_mixed_example()
  design <- compile_permanova_design(
    ~ treatment + (1 | cluster), example$data, dist(example$response)
  )
  fit <- permanova_dispersion(
    design, "treatment", type = "centroid", permutations = 19,
    seed = 73, keep_permutations = TRUE
  )
  cluster_means <- aggregate(
    fit$distances,
    list(cluster = example$data$cluster, treatment = example$data$treatment),
    mean
  )
  reference <- summary(aov(x ~ treatment, cluster_means))[[1]]

  expect_s3_class(fit, "permanova_dispersion")
  expect_equal(fit$engine, "balanced random intercept")
  expect_equal(
    fit$test$statistic,
    unname(reference["treatment", "F value"]),
    tolerance = 1e-10
  )
  expect_equal(fit$test$exchangeable_units, "cluster")
  expect_true(fit$test$symbolic_ems$exact_match)
})

test_that("nested dispersion uses its assigned hierarchy stratum", {
  example <- .dispersion_nested_example()
  design <- compile_permanova_design(
    ~ treatment + (1 | region/site),
    example$data, dist(example$response)
  )
  fit <- suppressWarnings(permanova_dispersion(
    design, "treatment", type = "centroid", permutations = "all",
    keep_permutations = TRUE
  ))
  site_means <- aggregate(
    fit$distances,
    list(region = example$data$region, site = example$data$site,
         treatment = example$data$treatment),
    mean
  )
  reference <- summary(aov(x ~ region + treatment, site_means))[[1]]

  expect_equal(fit$engine, "balanced nested hierarchy")
  expect_equal(fit$test$exchangeable_units, "region:site")
  expect_equal(fit$test$permutation_strata, "region")
  expect_equal(fit$test$attainable_permutations, 2 ^ 3)
  expect_equal(
    fit$test$statistic,
    unname(reference["treatment", "F value"]), tolerance = 1e-10
  )
})

test_that("crossed dispersion uses the treatment-by-block denominator", {
  example <- .dispersion_crossed_example()
  design <- compile_permanova_design(
    ~ treatment + (1 | block) + (1 | treatment:block),
    example$data, dist(example$response)
  )
  fit <- suppressWarnings(permanova_dispersion(
    design, "treatment", type = "centroid", permutations = "all",
    keep_permutations = TRUE
  ))
  reference_data <- data.frame(
    dispersion = fit$distances,
    treatment = example$data$treatment,
    block = example$data$block
  )
  table <- summary(aov(
    dispersion ~ treatment * block, reference_data
  ))[[1]]
  rows <- setNames(seq_len(nrow(table)), trimws(rownames(table)))
  reference <- table[rows[["treatment"]], "Mean Sq"] /
    table[rows[["treatment:block"]], "Mean Sq"]

  expect_equal(fit$engine, "balanced crossed randomized block")
  expect_equal(fit$test$exchangeable_units, "treatment:block")
  expect_equal(fit$test$permutation_strata, "block")
  expect_equal(fit$test$attainable_permutations, 2 ^ 3)
  expect_equal(fit$test$statistic, unname(reference), tolerance = 1e-10)
})

test_that("non-Euclidean dispersion retains the complete spectrum", {
  skip_if_not_installed("vegan")
  data(dune, dune.env, package = "vegan")
  complete <- dune.env$Management %in% c("NM", "SF")
  data <- droplevels(dune.env[complete, , drop = FALSE])
  distance <- as.dist(as.matrix(vegan::vegdist(dune))[complete, complete])
  design <- compile_permanova_design(~ Management, data, distance)
  fit <- permanova_dispersion(
    design, "Management", permutations = 19, seed = 75
  )

  expect_equal(fit$eigenvalues, design$spectral$eigenvalues, tolerance = 1e-10)
  expect_gt(fit$geometry$negative_count, 0)
  expect_true(is.finite(fit$test$p_value))
  expect_match(fit$test$warnings, "Negative principal-coordinate")
})

test_that("unsupported dispersion targets and structures fail explicitly", {
  dat <- data.frame(
    treatment = factor(rep(c("a", "b"), each = 4)),
    nuisance = factor(rep(c("x", "y"), 4)),
    continuous = seq_len(8), cluster = factor(rep(1:4, each = 2))
  )
  response <- matrix(seq_len(24), ncol = 3)

  multiple <- compile_permanova_design(
    ~ treatment + nuisance, dat, dist(response)
  )
  expect_error(
    permanova_dispersion(multiple, "treatment", permutations = 0),
    "exactly the one requested fixed term"
  )
  continuous <- compile_permanova_design(
    ~ continuous, dat, dist(response)
  )
  expect_error(
    permanova_dispersion(continuous, "continuous", permutations = 0),
    "categorical main-effect factor"
  )
  unbalanced <- dat[-1L, ]
  unbalanced_design <- compile_permanova_design(
    ~ treatment, unbalanced, dist(response[-1L, ])
  )
  expect_error(
    permanova_dispersion(unbalanced_design, "treatment", permutations = 0),
    "equal observations"
  )
  expect_error(
    permanova_dispersion(
      compile_permanova_design(
        ~ treatment + (1 | cluster), dat, dist(response)
      ),
      "nuisance", permutations = 0
    ),
    "one requested fixed term"
  )
})
