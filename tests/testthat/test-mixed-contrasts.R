.scalar_unit_contrast_f <- function(y, treatment, contrast, strata = NULL) {
  frame <- data.frame(treatment = treatment)
  if (!is.null(strata)) frame$strata <- strata
  X <- if (is.null(strata)) {
    model.matrix(~ treatment, frame)
  } else {
    model.matrix(~ strata + treatment, frame)
  }
  beta <- qr.solve(X, y)
  level_rows <- do.call(rbind, lapply(levels(treatment), function(level) {
    colMeans(X[treatment == level, , drop = FALSE])
  }))
  L <- as.numeric(contrast %*% level_rows)
  estimate <- sum(L * beta)
  numerator <- estimate ^ 2 / as.numeric(L %*% solve(crossprod(X), L))
  residual <- y - X %*% beta
  denominator <- sum(residual ^ 2) / (length(y) - qr(X)$rank)
  numerator / denominator
}

.three_level_contrast <- c(A = 1, B = -1, C = 0)

test_that("random-intercept contrasts use whole-unit centroids", {
  data <- expand.grid(replicate = factor(1:2), site = factor(1:9))
  data$treatment <- factor(rep(rep(c("A", "B", "C"), each = 3), each = 2))
  site_signal <- c(-1, 0, 1, 0, 1, 2, 1, 2, 4)
  data$y <- site_signal[as.integer(data$site)] + c(-.2, .2)
  design <- compile_permanova_design(
    ~ treatment + (1 | site), data, dist(data$y)
  )
  fit <- permanova_mixed_contrasts(
    design, .three_level_contrast, permutations = 0
  )
  unit <- aggregate(y ~ site + treatment, data, mean)
  reference <- .scalar_unit_contrast_f(
    unit$y, unit$treatment, .three_level_contrast
  )

  expect_equal(fit$table$pseudo_f, reference, tolerance = 1e-11)
  expect_match(fit$estimand, "random-intercept unit")
  expect_equal(fit$table$denominator_df, 6)
})

test_that("nested contrasts move intact top-level hierarchy units", {
  data <- expand.grid(
    replicate = factor(1:2), site = factor(1:2), region = factor(1:9)
  )
  data$treatment <- factor(rep(rep(c("A", "B", "C"), each = 3), each = 4))
  data$site_id <- interaction(data$region, data$site, drop = TRUE)
  region_signal <- c(-1, 0, 1, 0, 1, 2, 1, 2, 4)
  data$y <- region_signal[as.integer(data$region)] +
    .15 * as.integer(data$site) + c(-.1, .1)
  design <- compile_permanova_design(
    ~ treatment + (1 | region/site_id), data, dist(data$y)
  )
  fit <- permanova_mixed_contrasts(
    design, .three_level_contrast, permutations = 0
  )
  unit <- aggregate(y ~ region + treatment, data, mean)
  reference <- .scalar_unit_contrast_f(
    unit$y, unit$treatment, .three_level_contrast
  )

  expect_equal(fit$table$pseudo_f, reference, tolerance = 1e-10)
  expect_match(fit$design_type, "top-level nested")
  expect_equal(fit$permutation$units, "region")
})

test_that("crossed contrasts use the treatment-by-block denominator", {
  data <- expand.grid(
    replicate = factor(1:2), treatment = factor(c("A", "B", "C")),
    block = factor(1:5)
  )
  cell <- interaction(data$treatment, data$block, drop = TRUE)
  data$y <- rep(c(-1, 0, 1, .2, 1, 2, 1, 2, 4, -.5, 1, 3, 0, 2, 5), each = 2) +
    c(-.1, .1)
  design <- compile_permanova_design(
    ~ treatment + (1 | block) + (1 | treatment:block), data, dist(data$y)
  )
  fit <- permanova_mixed_contrasts(
    design, .three_level_contrast, permutations = 0
  )
  unit <- aggregate(y ~ treatment + block, data, mean)
  reference <- .scalar_unit_contrast_f(
    unit$y, unit$treatment, .three_level_contrast, unit$block
  )

  expect_equal(fit$table$pseudo_f, reference, tolerance = 1e-10)
  expect_equal(fit$denominator, "MS[treatment:block]")
  expect_equal(fit$table$denominator_df, 8)
  expect_equal(nlevels(cell), 15)
})

test_that("hybrid contrasts use treatment-by-region profiles", {
  data <- expand.grid(
    replicate = factor(1:2), site = factor(1:2),
    treatment = factor(c("A", "B", "C")), region = factor(1:5)
  )
  data$site_id <- interaction(data$region, data$site, drop = TRUE)
  profile <- interaction(data$treatment, data$region, drop = TRUE)
  profile_signal <- c(-1, 0, 1, .2, 1, 2, 1, 2, 4, -.5, 1, 3, 0, 2, 5)
  data$y <- profile_signal[profile] + .1 * as.integer(data$site) +
    c(-.05, .05)
  design <- compile_permanova_design(
    ~ treatment + (1 | region/site_id) + (1 | treatment:region) +
      (1 | treatment:region:site_id), data, dist(data$y)
  )
  fit <- permanova_mixed_contrasts(
    design, .three_level_contrast, permutations = 0
  )
  unit <- aggregate(y ~ treatment + region, data, mean)
  reference <- .scalar_unit_contrast_f(
    unit$y, unit$treatment, .three_level_contrast, unit$region
  )

  expect_equal(fit$table$pseudo_f, reference, tolerance = 1e-10)
  expect_equal(fit$denominator, "MS[treatment:region]")
  expect_match(fit$estimand, "profile centroids")
})

test_that("blocked mixed permutations equal direct response reconstruction", {
  data <- expand.grid(
    replicate = factor(1:2), treatment = factor(c("A", "B", "C")),
    block = factor(1:3)
  )
  data$y <- rep(c(-1, 0, 1, .2, 1, 2, 1, 2, 4), each = 2) + c(-.1, .1)
  design <- compile_permanova_design(
    ~ treatment + (1 | block) + (1 | treatment:block), data, dist(data$y)
  )
  fit <- permanova_mixed_contrasts(
    design, .three_level_contrast, permutations = "all",
    keep_permutations = TRUE
  )
  structure <- distmer:::.mixed_contrast_structure(design)
  units <- distmer:::.prepare_mixed_contrast_units(structure)
  plan <- distmer:::.mixed_contrast_plan(structure, "all", 100000)
  L <- as.numeric(.three_level_contrast %*% units$mean_rows)
  projector <- distmer:::.linear_contrast_projector(
    units$X, L, units$H_full, 1e-10
  )
  y <- vapply(structure$aggregated$blocks, function(index) mean(data$y[index]),
              numeric(1))
  residual <- diag(length(y)) - projector$H_reduced
  direct <- apply(plan$orders, 1, function(order) {
    y_permuted <- (projector$H_reduced + residual[order, , drop = FALSE]) %*% y
    (sum(y_permuted * (projector$H %*% y_permuted))) /
      (sum(y_permuted * (units$H_denominator %*% y_permuted)) /
         units$denominator_df)
  })

  expect_equal(fit$tests[[1]]$permutation_statistics, direct,
               tolerance = 1e-10)
  expect_true(fit$permutation$exact)
  expect_equal(fit$permutation$attainable, factorial(3)^3)
})

test_that("mixed pairwise adjustment and unsupported designs are explicit", {
  data <- expand.grid(replicate = factor(1:2), site = factor(1:9))
  data$treatment <- factor(rep(rep(c("A", "B", "C"), each = 3), each = 2))
  data$y <- as.integer(data$site) + c(-.1, .1)
  design <- compile_permanova_design(
    ~ treatment + (1 | site), data, dist(data$y)
  )
  pairwise <- permanova_mixed_pairwise(
    design, permutations = 19, seed = 4, adjust = "bonferroni"
  )

  expect_equal(nrow(pairwise$table), 3)
  expect_equal(pairwise$table$p_adjusted,
               p.adjust(pairwise$table$p_perm, "bonferroni"))

  fixed <- compile_permanova_design(~ treatment, data, dist(data$y))
  expect_error(permanova_mixed_contrasts(fixed, .three_level_contrast),
               "support one random")
  unequal <- data[-1, ]
  unequal_design <- compile_permanova_design(
    ~ treatment + (1 | site), unequal, dist(unequal$y)
  )
  unequal_fit <- permanova_mixed_contrasts(
    unequal_design, .three_level_contrast, permutations = 0
  )
  expect_match(unequal_fit$design_type, "unequally replicated")
})
