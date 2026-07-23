# P(MC) is the Monte-Carlo asymptotic p-value: it draws the pseudo-F from the
# eigenvalue-weighted chi-square distribution implied by the Gower matrix, so it
# stays smooth when the true permutation space is tiny. Its correctness anchor
# is the univariate Euclidean limit, where it must reproduce the exact classical
# F distribution -- for a one-way design and for a balanced mixed design.

test_that("P(MC) reproduces the classical F p-value for a univariate one-way", {
  set.seed(11)
  d <- data.frame(g = factor(rep(1:4, each = 12)))
  d$y <- rnorm(nrow(d)) + c(0, 0.6, 0, 0.3)[d$g]
  classical <- anova(aov(y ~ g, d))["g", "Pr(>F)"]
  fm <- permanova(dist(d$y) ~ g, d, permutations = 0,
                  mc_permutations = 200000, seed = 1)
  mc <- fm$table$`P(MC)`[fm$table$Source == "g"]
  expect_equal(mc, classical, tolerance = 0.01)
})

test_that("P(MC) reproduces the exact balanced mixed F-test over its stratum", {
  set.seed(101)
  d <- expand.grid(rep = 1:5, Site = factor(1:12))
  d$Loc <- factor(rep(c("A", "B", "C"), each = 20))
  d$y <- c(0, 0.9, 1.8)[d$Loc] + rnorm(12, 0, 1)[d$Site] + rnorm(nrow(d))
  fm <- permanova(dist(d$y) ~ Loc + (1 | Site), d, permutations = 0,
                  mc_permutations = 300000, seed = 1)
  row <- fm$table[fm$table$Source == "Loc", ]
  den_df <- fm$fit$table$denominator_df[fm$fit$table$term == "Loc"]
  classical <- pf(row$`Pseudo-F`, row$df, den_df, lower.tail = FALSE)
  expect_equal(row$`P(MC)`, classical, tolerance = 0.01)
})

test_that("a single non-zero eigenvalue collapses P(MC) exactly onto F", {
  # a univariate response has a rank-1 Gower matrix, so the internal draw is
  # lambda * chisq(df_num) over lambda * chisq(df_den): an exact F variate.
  set.seed(5)
  lambda <- 7.3
  f_obs <- 2.1
  p <- distmer:::.monte_carlo_p(lambda, df_num = 3, df_den = 24,
                                f_obs = f_obs, nmc = 500000)
  expect_equal(p, pf(f_obs, 3, 24, lower.tail = FALSE), tolerance = 0.005)
})

test_that("mc_permutations = 0 disables the column", {
  set.seed(2)
  d <- data.frame(g = factor(rep(1:3, each = 10)), y = rnorm(30))
  fm <- permanova(dist(d$y) ~ g, d, permutations = 0, mc_permutations = 0)
  expect_true(all(is.na(fm$table$`P(MC)`)))
})

test_that("P(MC) is reproducible from the seed and left NA for aggregated engines", {
  set.seed(3)
  d <- data.frame(g = factor(rep(1:3, each = 12)), y = rnorm(36))
  a <- permanova(dist(d$y) ~ g, d, permutations = 0, mc_permutations = 50000, seed = 42)
  b <- permanova(dist(d$y) ~ g, d, permutations = 0, mc_permutations = 50000, seed = 42)
  expect_identical(a$table$`P(MC)`, b$table$`P(MC)`)

  # an unequally replicated random-intercept design tests in a reduced
  # centroid space, so no full-G P(MC) is reported.
  set.seed(9)
  s <- data.frame(Site = factor(rep(1:8, times = c(3, 5, 4, 6, 3, 5, 4, 6))))
  s$Grp <- factor(c("A", "A", "A", "A", "B", "B", "B", "B")[as.integer(s$Site)])
  s$y <- rnorm(nrow(s)) + rnorm(8)[s$Site]
  fm <- permanova(dist(s$y) ~ Grp + (1 | Site), s,
                  permutations = 0, mc_permutations = 50000, seed = 1)
  expect_true(all(is.na(fm$table$`P(MC)`)))
})

test_that("mixed P(MC) is invariant to geometry outside the tested stratum", {
  set.seed(7301)
  ng <- 12L; m <- 5L; q <- 5L
  d <- expand.grid(rep = seq_len(m), Site = factor(seq_len(ng)))
  d$Loc <- factor(rep(rep(c("A", "B", "C"), each = 4L), each = m))
  centroids <- matrix(rnorm(ng * q), ng, q)
  deviations <- matrix(rnorm(nrow(d) * q), nrow(d), q)
  deviations <- do.call(rbind, lapply(
    split(seq_len(nrow(d)), d$Site),
    function(index) scale(deviations[index, , drop = FALSE],
                          center = TRUE, scale = FALSE)
  ))
  response <- function(scale) {
    centroids[as.integer(d$Site), , drop = FALSE] + scale * deviations
  }
  small <- permanova(
    dist(response(0.01)) ~ Loc + (1 | Site), d,
    permutations = 0, mc_permutations = 50000, seed = 91
  )
  large <- permanova(
    dist(response(100)) ~ Loc + (1 | Site), d,
    permutations = 0, mc_permutations = 50000, seed = 91
  )
  expect_equal(small$table$`Pseudo-F`[1L], large$table$`Pseudo-F`[1L],
               tolerance = 1e-8)
  expect_identical(small$table$`P(MC)`[1L], large$table$`P(MC)`[1L])
})
