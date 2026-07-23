# A continuous random slope aggregates to one within-group slope vector per
# group. For a balanced Euclidean design this is exactly the two-stage slope
# test: fit a slope per group, one-sample test on the group slopes, with
# reference F(1, ng - 1). That exact identity is the correctness anchor.

.slope_example <- function(seed = 21, ng = 14, m = 8) {
  set.seed(seed)
  d <- data.frame(g = factor(rep(seq_len(ng), each = m)))
  d$x <- rep(seq(-2, 2, length.out = m), ng)
  b0 <- rnorm(ng, 0, 1.5)[d$g]; b1 <- rnorm(ng, 0, 0.7)[d$g]
  d$y <- 1 + 0.5 * d$x + b0 + b1 * d$x + rnorm(nrow(d), 0, 1)
  d
}

test_that("continuous random slope equals the exact two-stage slope test", {
  d <- .slope_example()
  ng <- nlevels(d$g)
  slopes <- vapply(split(d, d$g), function(s) coef(lm(y ~ x, s))[[2L]], numeric(1))
  f_two <- ng * mean(slopes)^2 / var(slopes)
  p_two <- pf(f_two, 1, ng - 1, lower.tail = FALSE)

  fm <- permanova_random_slope_continuous(
    compile_permanova_design(~ x + (1 + x | g), d, dist(d$y)), permutations = 0
  )
  expect_equal(fm$test$statistic, f_two, tolerance = 1e-8)
  expect_equal(fm$test$p_value_parametric, p_two, tolerance = 1e-8)
  expect_equal(fm$test$denominator$df, ng - 1L)
})

test_that("the test is invariant to within-group centring of the covariate", {
  d <- .slope_example()
  d$xw <- decompose_covariate(d, "x", by = "g")$x_within
  raw <- permanova_random_slope_continuous(
    compile_permanova_design(~ x + (1 + x | g), d, dist(d$y)), permutations = 0
  )
  cen <- permanova_random_slope_continuous(
    compile_permanova_design(~ xw + (1 + xw | g), d, dist(d$y)), permutations = 0
  )
  expect_equal(raw$test$statistic, cen$test$statistic, tolerance = 1e-8)
})

test_that("permanova() routes a numeric random slope to the continuous engine", {
  d <- .slope_example()
  fm <- permanova(dist(d$y) ~ x + (1 + x | g), d, permutations = 199, seed = 1)
  expect_identical(fm$engine, "continuous random slope")
  # the aggregated slope space does not partition n - 1, so no Total / P(MC)
  expect_false(any(fm$table$Source == "Total"))
  expect_true(all(is.na(fm$table$`P(MC)`)))
})

test_that("whole-group sign flips give a permutation p-value near the exact one", {
  d <- .slope_example()
  fm <- permanova_random_slope_continuous(
    compile_permanova_design(~ x + (1 + x | g), d, dist(d$y)),
    permutations = 4999, seed = 3
  )
  expect_lt(abs(fm$test$p_value - fm$test$p_value_parametric), 0.03)
})

test_that("a non-Euclidean dissimilarity suppresses the parametric p-value only", {
  set.seed(5)
  ng <- 12; m <- 6
  d <- data.frame(g = factor(rep(seq_len(ng), each = m)))
  d$x <- rnorm(nrow(d))
  Y <- matrix(rpois(nrow(d) * 5, 4) + 1, nrow(d), 5)
  Y <- Y + round(outer(rnorm(ng)[d$g] * d$x, 1:5))
  Y[Y < 0] <- 0
  fm <- permanova_random_slope_continuous(
    compile_permanova_design(~ x + (1 + x | g), d, vegan::vegdist(Y, "bray")),
    permutations = 999, seed = 1
  )
  expect_true(is.na(fm$test$p_value_parametric))
  expect_false(is.na(fm$test$p_value))
})

test_that("the engine rejects too few groups and inestimable within-group slopes", {
  d <- .slope_example(ng = 2)
  expect_error(
    permanova_random_slope_continuous(
      compile_permanova_design(~ x + (1 + x | g), d, dist(d$y)), permutations = 0
    ), "three groups"
  )
  d2 <- .slope_example(ng = 5, m = 6)
  d2$x[d2$g == "1"] <- 4  # no within-group variation in group 1
  expect_error(
    permanova_random_slope_continuous(
      compile_permanova_design(~ x + (1 + x | g), d2, dist(d2$y)), permutations = 0
    ), "not estimable"
  )
})

test_that("parametric slope inference is suppressed for unequal information", {
  set.seed(22)
  sizes <- c(3, 4, 5, 7, 10, 15)
  d <- data.frame(g = factor(rep(seq_along(sizes), sizes)))
  d$x <- unlist(Map(
    function(n, scale) seq(-scale, scale, length.out = n),
    sizes, c(0.1, 0.2, 0.5, 1, 2, 4)
  ))
  slope <- rnorm(length(sizes), 0, 0.8)
  d$y <- slope[d$g] * d$x + rnorm(nrow(d))
  fm <- permanova_random_slope_continuous(
    compile_permanova_design(~ x + (1 + x | g), d, dist(d$y)),
    permutations = 0
  )
  expect_true(is.na(fm$test$p_value_parametric))
  expect_true(any(grepl("SSx differs", fm$test$warnings, fixed = TRUE)))

  slopes <- vapply(split(seq_len(nrow(d)), d$g), function(index) {
    coef(lm(d$y[index] ~ d$x[index]))[[2L]]
  }, numeric(1))
  ssx <- vapply(split(d$x, d$g), function(x) sum((x - mean(x))^2), numeric(1))
  residual <- fm$variance_components$unrestricted[["residual"]]
  reference <- var(slopes) - residual * mean(1 / ssx)
  expect_equal(fm$variance_components$unrestricted[["slope"]], reference,
               tolerance = 1e-8)
})

test_that("multivariate Euclidean slopes do not receive a univariate F reference", {
  d <- .slope_example()
  set.seed(12)
  Y <- cbind(d$y, 0.4 * d$x + rnorm(nrow(d)))
  fm <- permanova_random_slope_continuous(
    compile_permanova_design(~ x + (1 + x | g), d, dist(Y)),
    permutations = 0
  )
  expect_true(is.na(fm$test$p_value_parametric))
  expect_true(any(grepl("not rank one", fm$test$warnings, fixed = TRUE)))
})
