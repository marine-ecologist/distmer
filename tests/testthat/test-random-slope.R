.random_slope_example <- function(subjects = 6L, replicates = 2L) {
  data <- expand.grid(
    replicate = seq_len(replicates),
    condition = factor(c("before", "after"), levels = c("before", "after")),
    subject = factor(seq_len(subjects))
  )
  intercept <- seq(-0.8, 0.7, length.out = subjects)
  slope <- c(-0.4, 0.1, 0.3, 0.8, 1.0, 1.4)[seq_len(subjects)]
  data$y <- intercept[data$subject] +
    slope[data$subject] * (data$condition == "after") +
    rep(c(-0.12, 0.12), length.out = nrow(data))
  data
}

.random_slope_formula <- function() {
  ~ condition + (1 + condition | subject)
}

.subject_coefficients <- function(response, data) {
  cell <- aggregate(
    response,
    list(subject = data$subject, condition = data$condition), mean
  )
  wide <- reshape(
    cell, idvar = "subject", timevar = "condition", direction = "wide"
  )
  cbind(
    midpoint = (wide$x.before + wide$x.after) / 2,
    difference = wide$x.after - wide$x.before
  )
}

test_that("random-slope F equals paired-difference t squared", {
  data <- .random_slope_example()
  design <- compile_permanova_design(
    .random_slope_formula(), data, dist(data$y)
  )
  fit <- suppressWarnings(permanova_random_slope(
    design, permutations = 0, truncate = TRUE
  ))
  coefficients <- .subject_coefficients(data$y, data)
  reference <- mean(coefficients[, "difference"])^2 /
    (stats::var(coefficients[, "difference"]) / nrow(coefficients))

  expect_s3_class(fit, "permanova_random_slope")
  expect_equal(fit$test$statistic, reference, tolerance = 1e-10)
  expect_equal(fit$test$statistic,
               unname(stats::t.test(coefficients[, "difference"])$statistic)^2,
               tolerance = 1e-10)
  expect_equal(fit$test$denominator$df, 5)
  expect_equal(fit$structure$replicates_per_cell, 2)
  expect_true(all(eigen(fit$covariance$truncated, symmetric = TRUE)$values >=
                    -1e-12))
})

test_that("random-slope covariance equals scalar method of moments", {
  data <- .random_slope_example()
  design <- compile_permanova_design(
    .random_slope_formula(), data, dist(data$y)
  )
  estimate <- estimate_random_slope_covariance(design)
  coefficients <- .subject_coefficients(data$y, data)
  cell <- interaction(data$subject, data$condition)
  residual <- data$y - ave(data$y, cell)
  residual_ms <- sum(residual^2) / (nrow(data) - nlevels(cell))
  empirical <- stats::cov(coefficients)
  reference <- empirical
  reference[1, 1] <- reference[1, 1] - residual_ms / 4
  reference[2, 2] <- reference[2, 2] - residual_ms
  dimnames(reference) <- dimnames(estimate$unrestricted)

  expect_equal(estimate$residual_trace_variance, residual_ms,
               tolerance = 1e-10)
  expect_equal(estimate$unrestricted, reference, tolerance = 1e-10)
  expect_match(estimate$parameterization, "after - before")
})

test_that("exact random-slope sign flips match direct difference vectors", {
  data <- .random_slope_example(subjects = 5L)
  design <- compile_permanova_design(
    .random_slope_formula(), data, dist(data$y)
  )
  fit <- suppressWarnings(permanova_random_slope(
    design, permutations = "all", keep_permutations = TRUE
  ))
  difference <- .subject_coefficients(data$y, data)[, "difference"]
  plan <- distmer:::.sign_flip_plan(5L, "all", 100000L)
  direct <- apply(plan$signs, 1L, function(sign) {
    value <- sign * difference
    mean(value)^2 / (stats::var(value) / length(value))
  })

  expect_true(fit$test$exact_permutation)
  expect_equal(fit$test$attainable_permutations, 2^5)
  expect_equal(fit$test$permutation_statistics, direct, tolerance = 1e-10)
  expect_equal(fit$test$p_value, mean(direct >= fit$test$statistic - 1e-10))
  expect_match(fit$test$exchangeability_assumption, "central symmetry")
})

test_that("random-slope Monte Carlo sign flips are seed reproducible", {
  data <- .random_slope_example()
  design <- compile_permanova_design(
    .random_slope_formula(), data, dist(data$y)
  )
  first <- permanova_random_slope(
    design, permutations = 9, seed = 27, keep_permutations = TRUE
  )
  second <- permanova_random_slope(
    design, permutations = 9, seed = 27, keep_permutations = TRUE
  )

  expect_false(first$test$exact_permutation)
  expect_equal(first$test$permutation_statistics,
               second$test$permutation_statistics)
  expect_equal(first$test$p_value, second$test$p_value)
})

test_that("unsupported random-slope envelopes reject explicitly", {
  data <- .random_slope_example()

  three <- data
  three$condition <- factor(rep(c("before", "during", "after"),
                                length.out = nrow(three)))
  three_design <- compile_permanova_design(
    ~ condition + (1 + condition | subject), three, dist(three$y)
  )
  expect_error(
    permanova_random_slope(three_design, permutations = 0), "two-level"
  )

  missing <- data[!(data$subject == "1" & data$condition == "after"), ]
  missing_design <- compile_permanova_design(
    .random_slope_formula(), missing, dist(missing$y)
  )
  expect_error(
    permanova_random_slope(missing_design, permutations = 0),
    "Every subject must contain both"
  )

  unequal <- data[-1L, ]
  unequal_design <- compile_permanova_design(
    .random_slope_formula(), unequal, dist(unequal$y)
  )
  expect_error(
    permanova_random_slope(unequal_design, permutations = 0),
    "equal replication"
  )

  unreplicated <- .random_slope_example(replicates = 1L)
  unreplicated_design <- compile_permanova_design(
    .random_slope_formula(), unreplicated, dist(unreplicated$y)
  )
  expect_error(
    permanova_random_slope(unreplicated_design, permutations = 0),
    "at least two observations"
  )

  extra <- compile_permanova_design(
    ~ condition + replicate + (1 + condition | subject), data, dist(data$y)
  )
  expect_error(
    permanova_random_slope(extra, permutations = 0), "sole fixed main effect"
  )
})

test_that("random-slope zero denominators reject and semimetrics stay empirical", {
  data <- .random_slope_example()
  constant_difference <- as.numeric(data$subject) +
    0.5 * (data$condition == "after")
  zero <- compile_permanova_design(
    .random_slope_formula(), data, dist(constant_difference)
  )
  expect_error(
    permanova_random_slope(zero, permutations = 0), "denominator is numerically zero"
  )

  distance <- as.matrix(dist(cbind(data$y, seq_len(nrow(data)) / 10)))
  distance[1, nrow(distance)] <- distance[nrow(distance), 1] <-
    5 * max(distance)
  semimetric <- compile_permanova_design(
    .random_slope_formula(), data, distance
  )
  fit <- permanova_random_slope(semimetric, permutations = 9, seed = 8)
  expect_false(semimetric$spectral$euclidean)
  expect_true(is.finite(fit$test$p_value))
  expect_true(is.na(fit$test$p_value_parametric))
  expect_match(paste(fit$test$warnings, collapse = " "), "suppressed")
})
