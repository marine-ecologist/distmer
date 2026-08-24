# Behaviour at extreme statistic magnitudes.
#
# Every engine is otherwise tested on well-conditioned data. The absolute
# comparison tolerance that produced p = 0 survived precisely because nothing
# exercised a large pseudo-F, so these checks pin the invariants that hold
# whatever the conditioning:
#   * a permutation p-value lies in (0, 1];
#   * under exact enumeration p >= 1/N, because the identity permutation is a
#     member of the set and always ties with the observed statistic;
#   * pseudo-F is scale free, so rescaling the response must not move a p-value.

.conditioning_balanced <- function(separation, seed = 7L) {
  d <- expand.grid(rep = 1:3, unit = factor(1:6))
  d$tr <- factor(ifelse(as.integer(d$unit) <= 3L, "control", "impact"))
  set.seed(seed)
  d$y <- separation * as.integer(d$tr) +
    stats::rnorm(nrow(d), 0, separation * 1e-4)
  d
}

.conditioning_unbalanced <- function(separation, seed = 7L) {
  d <- do.call(rbind, lapply(1:6, function(u) {
    data.frame(unit = paste0("u", u), rep = seq_len(if (u %% 2L) 2L else 3L))
  }))
  d$unit <- factor(d$unit)
  d$tr <- factor(ifelse(as.integer(d$unit) <= 3L, "control", "impact"))
  set.seed(seed)
  d$y <- separation * as.integer(d$tr) +
    stats::rnorm(nrow(d), 0, separation * 1e-4)
  d
}

test_that("exact permutation p-values stay positive at extreme separation", {
  for (separation in c(1, 1e3)) {
    balanced <- .conditioning_balanced(separation)
    unbalanced <- .conditioning_unbalanced(separation)

    fits <- list(
      mixed = suppressWarnings(permanova_mixed(
        compile_permanova_design(~ tr + (1 | unit), balanced, dist(balanced$y)),
        permutations = "all")),
      unequal = permanova_unequal_random_intercept(
        compile_permanova_design(~ tr + (1 | unit), unbalanced, dist(unbalanced$y)),
        permutations = "all"),
      splitplot = permanova_unbalanced_splitplot(
        compile_permanova_design(~ tr + (1 | unit), unbalanced, dist(unbalanced$y)),
        permutations = "all")
    )
    for (name in names(fits)) {
      fit <- fits[[name]]
      p <- fit$table$p_perm[[1L]]
      used <- fit$permutation$used
      expect_gt(p, 0)
      expect_lte(p, 1)
      # the identity permutation guarantees this floor
      expect_gte(p, 1 / used - 1e-12)
    }
  }
})

test_that("the exact p-value floor is exactly 1/N when ties are lost to rounding", {
  # a separation severe enough that the observed statistic and its own permuted
  # counterpart no longer agree to any fixed relative slack
  d <- .conditioning_unbalanced(1e3)
  fit <- permanova_unequal_random_intercept(
    compile_permanova_design(~ tr + (1 | unit), d, dist(d$y)),
    permutations = "all"
  )
  expect_gt(fit$table$pseudo_f[[1L]], 1e7)
  expect_equal(fit$table$p_perm[[1L]], 1 / fit$permutation$used)
})

test_that("the floor applies to the shared helper and only when exact", {
  statistics <- rep(0.1, 20L)          # nothing reaches the observed value
  exact <- .permutation_inference(statistics, observed = 99, exact = TRUE)
  expect_equal(exact$p_value, 1 / 20)

  # Monte Carlo already carries its own +1 and must not be floored twice
  sampled <- .permutation_inference(statistics, observed = 99, exact = FALSE)
  expect_equal(sampled$p_value, 1 / 21)
})

test_that("p-values do not move when the response is rescaled", {
  ordinary <- data.frame(g = factor(rep(c("a", "b", "c"), each = 4L)))
  set.seed(11)
  ordinary$signal <- as.numeric(ordinary$g) + stats::rnorm(nrow(ordinary), 0, 0.3)
  # ordinary noise, so rescaling does not drive the residual mean square under
  # the (absolute) numerically-zero guard the engines still apply
  balanced <- expand.grid(rep = 1:3, unit = factor(1:6))
  balanced$tr <- factor(ifelse(as.integer(balanced$unit) <= 3L, "control", "impact"))
  set.seed(7)
  balanced$y <- 2 * as.integer(balanced$tr) +
    stats::rnorm(nlevels(balanced$unit))[balanced$unit] +
    stats::rnorm(nrow(balanced), 0, 0.5)

  engines <- list(
    oneway = function(y) permanova_oneway(dist(y), ordinary$g,
                                          permutations = 199L, seed = 1L)$table$p_perm[[1L]],
    unbalanced = function(y) permanova_unbalanced(
      compile_permanova_design(~ g, ordinary, dist(y)),
      permutations = 199L, seed = 1L)$table$p_perm[[1L]],
    maxT = function(y) permanova_maxT_pairwise(
      compile_permanova_design(~ g, ordinary, dist(y)),
      permutations = 199L, seed = 1L)$table$p_perm[[1L]],
    dispersion = function(y) permanova_dispersion(
      compile_permanova_design(~ g, ordinary, dist(y)),
      permutations = 199L, seed = 1L)$table$p_perm[[1L]]
  )
  for (name in names(engines)) {
    values <- vapply(c(1e-4, 1, 1e4),
                     function(s) engines[[name]](ordinary$signal * s), numeric(1))
    expect_equal(values[[1L]], values[[2L]], info = name)
    expect_equal(values[[2L]], values[[3L]], info = name)
  }

  mixed <- vapply(c(1e-4, 1, 1e4), function(s) {
    suppressWarnings(permanova_mixed(
      compile_permanova_design(~ tr + (1 | unit), balanced, dist(balanced$y * s)),
      permutations = "all"))$table$p_perm[[1L]]
  }, numeric(1))
  expect_equal(mixed[[1L]], mixed[[2L]])
  expect_equal(mixed[[2L]], mixed[[3L]])
})

test_that("a residual negligible against total inertia is refused", {
  # Separation so severe that the residual is ~1e-12 of the total inertia: the
  # ratio being tested is no longer meaningful, and the engine must refuse
  # rather than report a pseudo-F built from rounding noise.
  d <- expand.grid(rep = 1:3, unit = factor(1:6))
  d$tr <- factor(ifelse(as.integer(d$unit) <= 3L, "control", "impact"))
  set.seed(7)
  d$y <- 1e3 * as.integer(d$tr) + stats::rnorm(nrow(d), 0, 1e-3)
  expect_error(
    suppressWarnings(permanova_mixed(
      compile_permanova_design(~ tr + (1 | unit), d, dist(d$y)),
      permutations = "all")),
    "numerically zero"
  )
})

test_that("the numerically-zero guard is scale free", {
  d <- expand.grid(rep = 1:3, unit = factor(1:6))
  d$tr <- factor(ifelse(as.integer(d$unit) <= 3L, "control", "impact"))
  set.seed(7)
  d$y <- 2 * as.integer(d$tr) +
    stats::rnorm(nlevels(d$unit))[d$unit] + stats::rnorm(nrow(d), 0, 0.5)
  # the same analysis must be accepted or refused identically at any scale;
  # before 0.9.4 the absolute guard rejected the small-scale copy
  for (scale in c(1e-6, 1e-3, 1, 1e3, 1e6)) {
    expect_no_error(suppressWarnings(permanova_mixed(
      compile_permanova_design(~ tr + (1 | unit), d, dist(d$y * scale)),
      permutations = 0)))
  }
})

test_that("a retained non-Euclidean spectrum is not mistaken for no variation", {
  # the aggregated Gower matrix here has a negative trace; that is a supported
  # result, not an absence of variation, so the guard must use magnitudes
  G <- diag(c(3.7, 1.0, 0, 0, 0, -50.9))
  expect_false(.mean_square_is_zero(-10.1, G))
  expect_true(.mean_square_is_zero(1e-14, G))
})
