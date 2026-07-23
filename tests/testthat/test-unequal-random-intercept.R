.unequal_cluster_example <- function() {
  sizes <- c(2L, 2L, 3L, 3L, 4L, 4L)
  dat <- data.frame(
    site = factor(rep(seq_along(sizes), sizes)),
    treatment = factor(rep(rep(c("control", "impact"), 3L), sizes))
  )
  site_effect <- c(-0.8, 0.4, -0.2, 0.7, -0.5, 1.0)
  dat$y <- 0.6 * (dat$treatment == "impact") +
    site_effect[dat$site] + ave(seq_len(nrow(dat)), dat$site, FUN = function(x) {
      seq(-0.15, 0.15, length.out = length(x))
    })
  dat
}

test_that("unequal random-intercept F equals ANOVA on cluster centroids", {
  dat <- .unequal_cluster_example()
  design <- compile_permanova_design(
    ~ treatment + (1 | site), dat, dist(dat$y)
  )
  fit <- suppressWarnings(permanova_unequal_random_intercept(
    design, permutations = "all", keep_permutations = TRUE
  ))
  means <- aggregate(y ~ site + treatment, dat, mean)
  reference <- summary(aov(y ~ treatment, means))[[1]][1L, "F value"]

  expect_equal(fit$test$statistic, unname(reference), tolerance = 1e-10)
  expect_equal(fit$structure$unit_sizes, c(2L, 2L, 3L, 3L, 4L, 4L))
  expect_equal(fit$permutation$attainable, 8)
  expect_equal(fit$permutation$used, 8)
  expect_equal(fit$estimand, "equally weighted random-intercept unit centroids")
})

test_that("size-stratified residual permutations match a direct scalar fit", {
  dat <- .unequal_cluster_example()
  design <- compile_permanova_design(
    ~ treatment + (1 | site), dat, dist(dat$y)
  )
  fit <- suppressWarnings(permanova_unequal_random_intercept(
    design, permutations = "all", keep_permutations = TRUE
  ))
  structure <- .validate_unequal_random_intercept(design)
  plan <- .stratified_permutation_plan(structure$size_strata, "all", 100L)
  means <- as.numeric(aggregate(dat$y, list(dat$site), mean)$x)
  reduced_fit <- mean(means)
  residual <- means - reduced_fit
  direct <- vapply(seq_len(nrow(plan$orders)), function(row) {
    permuted <- reduced_fit + residual[plan$orders[row, ]]
    table <- summary(aov(permuted ~ structure$treatment))[[1]]
    unname(table[1L, "F value"])
  }, numeric(1))
  expect_equal(fit$test$permutation_statistics, direct, tolerance = 1e-10)
})

test_that("Phase 11 rejects unsafe unit definitions", {
  dat <- .unequal_cluster_example()
  varying <- dat
  varying$treatment[1L] <- "impact"
  varying_design <- compile_permanova_design(
    ~ treatment + (1 | site), varying, dist(varying$y)
  )
  expect_error(
    permanova_unequal_random_intercept(varying_design, permutations = 0),
    "constant within"
  )

  confounded <- dat
  confounded$treatment <- factor(ifelse(
    ave(seq_len(nrow(confounded)), confounded$site, FUN = length) <= 2,
    "control", "impact"
  ))
  confounded_design <- compile_permanova_design(
    ~ treatment + (1 | site), confounded, dist(confounded$y)
  )
  expect_error(
    permanova_unequal_random_intercept(confounded_design, permutations = 0),
    "do not span every treatment contrast"
  )

  dat$subsite <- factor(ave(
    seq_len(nrow(dat)), dat$site, FUN = function(x) rep(1:2, length.out = length(x))
  ))
  nested <- compile_permanova_design(
    ~ treatment + (1 | site/subsite), dat, dist(dat$y)
  )
  expect_error(
    permanova_unequal_random_intercept(nested, permutations = 0),
    "exactly one"
  )
})

test_that("Phase 11 rejects exact size-treatment confounding despite centering dust", {
  sizes <- c(2L, 2L, 2L, 3L, 3L, 3L)
  data <- data.frame(
    site = factor(rep(seq_along(sizes), sizes)),
    treatment = factor(rep(
      c("control", "control", "control", "impact", "impact", "impact"),
      sizes
    ))
  )
  data$y <- as.integer(data$site) +
    ave(seq_len(nrow(data)), data$site, FUN = function(index) {
      seq(-0.1, 0.1, length.out = length(index))
    })
  design <- compile_permanova_design(
    ~ treatment + (1 | site), data, dist(data$y)
  )

  X <- model.matrix(~ factor(c("control", "control", "control",
                               "impact", "impact", "impact")))
  rank <- distmer:::.stratified_treatment_rank(
    X[, -1L, drop = FALSE], factor(sizes), 1e-10
  )
  expect_equal(as.integer(rank), 0L)
  expect_lt(attr(rank, "projected_singular_values")[[1L]],
            attr(rank, "absolute_floor"))
  expect_error(
    permanova_unequal_random_intercept(design, permutations = 0),
    "do not span every treatment contrast"
  )
})

test_that("Phase 11 retains non-Euclidean empirical inference", {
  dat <- .unequal_cluster_example()
  D <- as.matrix(dist(cbind(dat$y, seq_len(nrow(dat)) / 10)))
  D[1L, 2L] <- D[2L, 1L] <- 5 * max(D)
  design <- compile_permanova_design(
    ~ treatment + (1 | site), dat, D
  )
  fit <- suppressWarnings(permanova_unequal_random_intercept(
    design, permutations = "all"
  ))
  expect_false(design$spectral$euclidean)
  expect_true(is.finite(fit$test$p_value))
  expect_true(any(grepl("non-Euclidean", fit$test$warnings)))
})
