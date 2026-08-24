# Extremeness comparison and permutation resolution.
#
# Two related failure modes live here. A restricted scheme can attain many
# arrangements but few distinct statistics, putting a floor on the p-value; and
# comparing statistics with an absolute tolerance breaks down once the statistic
# is large, because floating-point noise then exceeds the tolerance.

test_that("the p-value matches the classical idiom on ordinary scales", {
  set.seed(2)
  for (i in seq_len(200L)) {
    statistics <- round(stats::rnorm(50L), 2L)
    observed <- sample(statistics, 1L)
    exact <- sample(c(TRUE, FALSE), 1L)
    extreme <- sum(statistics >= observed - 1e-10, na.rm = TRUE)
    classical <- if (exact) {
      extreme / length(statistics)
    } else (extreme + 1) / (length(statistics) + 1)
    expect_equal(.permutation_inference(statistics, observed, exact)$p_value,
                 classical)
  }
})

test_that("a large statistic still counts its own permutation as extreme", {
  # The identity permutation reproduces the observed statistic up to rounding.
  # At a pseudo-F of ~1e4 that rounding is ~1e-9, which an absolute 1e-10
  # tolerance cannot absorb -- the exact test then reported p = 0.
  observed <- 9254.4400000045
  statistics <- c(9254.4400000019, rep(0.5, 23L))

  expect_gt(observed - max(statistics), 1e-10)          # the gap that broke it
  result <- .permutation_inference(statistics, observed, exact = TRUE)
  expect_equal(result$p_value, 1 / 24)
  expect_gt(result$p_value, 0)
})

test_that("an exact permutation p-value is never zero", {
  set.seed(5)
  d <- data.frame(unit = factor(rep(paste0("u", 1:4), each = 3L)))
  d$tr <- factor(ifelse(as.integer(d$unit) <= 2L, "control", "impact"))
  # enormous separation: every permuted statistic is tiny except the identity
  d$y <- c(1, 1.1, 0.9, 1.2, 1.0, 1.1, 9, 9.2, 8.8, 9.1, 9.3, 9.0)

  fit <- permanova_unequal_random_intercept(
    compile_permanova_design(~ tr + (1 | unit), d, dist(d$y)),
    permutations = "all"
  )
  expect_gt(fit$table$pseudo_f[[1L]], 1000)
  expect_gt(fit$table$p_perm[[1L]], 0)
  expect_equal(fit$table$p_perm[[1L]], 1 / 3, tolerance = 1e-8)
})

test_that("resolution is reported and a floored design is flagged", {
  # three distinct values, eight tied at the maximum out of 24
  statistics <- rep(c(0.1, 1.0, 5.0), times = c(8L, 8L, 8L))
  result <- .permutation_inference(statistics, observed = 5.0, exact = TRUE)

  expect_equal(result$distinct_statistics, 3L)
  expect_equal(result$minimum_attainable_p, 1 / 3)
  expect_match(result$resolution_note, "cannot reach 0.05")

  # a well-resolved distribution carries no note
  set.seed(1)
  healthy <- .permutation_inference(stats::rnorm(999L), observed = 4,
                                    exact = FALSE)
  expect_length(healthy$resolution_note, 0L)
  expect_gt(healthy$distinct_statistics, 900L)
})

test_that("the resolution warning reaches the fitted object", {
  set.seed(5)
  d <- data.frame(unit = factor(rep(paste0("u", 1:4), each = 3L)))
  d$tr <- factor(ifelse(as.integer(d$unit) <= 2L, "control", "impact"))
  d$y <- c(1, 1.1, 0.9, 1.2, 1.0, 1.1, 9, 9.2, 8.8, 9.1, 9.3, 9.0)
  fit <- permanova_unequal_random_intercept(
    compile_permanova_design(~ tr + (1 | unit), d, dist(d$y)),
    permutations = "all"
  )
  expect_true(any(grepl("distinct statistic", fit$test$warnings)))
})

test_that("an empty permutation set yields no p-value and no note", {
  result <- .permutation_inference(numeric(0), observed = 1, exact = TRUE)
  expect_true(is.na(result$p_value))
  expect_equal(result$distinct_statistics, 0L)
  expect_length(result$resolution_note, 0L)
})

test_that("distinct-value counting is tolerance aware", {
  expect_equal(.distinct_values(c(1, 1 + 1e-14, 2, 3), 1e-10), 3L)
  expect_equal(.distinct_values(c(1, 1.5, 2), 1e-10), 3L)
  expect_equal(.distinct_values(numeric(0), 1e-10), 0L)
  # counting must scale: values near 1e4 separated by 1e-9 are one value
  expect_equal(.distinct_values(c(9254.44, 9254.44 + 1e-9), 1e-10), 1L)
})

test_that("the reporting level is settable and validated", {
  statistics <- rep(c(0.1, 1.0, 5.0), times = c(8L, 8L, 8L))
  old <- options(distmer.resolution_alpha = NULL)
  on.exit(options(old), add = TRUE)

  # default 0.05: a design floored at 1/3 is flagged
  expect_length(
    .permutation_inference(statistics, 5.0, exact = TRUE)$resolution_note, 1L)

  # raise the level above the floor and the note is no longer warranted
  options(distmer.resolution_alpha = 0.5)
  quiet <- .permutation_inference(statistics, 5.0, exact = TRUE)
  expect_length(quiet$resolution_note, 0L)
  # the underlying numbers are recorded either way
  expect_equal(quiet$distinct_statistics, 3L)
  expect_equal(quiet$minimum_attainable_p, 1 / 3)

  options(distmer.resolution_alpha = 0.25)
  expect_match(.permutation_inference(statistics, 5.0, exact = TRUE)$resolution_note,
               "cannot reach 0.25")

  # the level never changes the p-value
  options(distmer.resolution_alpha = 0.5)
  expect_equal(.permutation_inference(statistics, 5.0, exact = TRUE)$p_value,
               1 / 3)

  for (bad in list(0, 1, -0.1, c(0.05, 0.1), NA_real_)) {
    options(distmer.resolution_alpha = bad)
    expect_error(.permutation_inference(statistics, 5.0, exact = TRUE),
                 "strictly between")
  }
})
