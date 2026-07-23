test_that("exact unique assignments agree with exhaustive full-engine permutations", {
  audit <- audit_unequal_random_intercept_calibration(
    batches = 2, simulations_per_batch = 10, permutations = 19,
    engine_checks = 3, seed = 20260719, keep_results = TRUE
  )
  expect_s3_class(audit, "permanova_random_intercept_calibration_audit")
  expect_lt(audit$max_engine_exact_f_difference, 1e-10)
  expect_lt(audit$max_engine_exact_p_difference, 1e-12)
  expect_equal(diff(audit$seed_design$batch_seeds), 1000000L)
  expect_equal(
    unique(audit$results$exact_assignments[
      audit$results$scenario == "gaussian_balanced"
    ]),
    216
  )
  expect_equal(
    unique(audit$results$exact_assignments[
      audit$results$scenario == "gaussian_unequal_allocation"
    ]),
    96
  )
})

test_that("intercept-only Freedman-Lane equals raw stratified permutation", {
  strata <- factor(rep(1:3, each = 4))
  plan <- .stratified_permutation_plan(strata, 5, 1000L)
  H0 <- matrix(1 / 12, 12, 12)
  R0 <- diag(12) - H0
  errors <- apply(plan$orders, 1L, function(index) {
    P <- diag(12)[index, , drop = FALSE]
    max(abs(H0 + R0[index, , drop = FALSE] - P))
  })
  expect_lt(max(errors), 1e-12)
})

test_that("Phase 11 distinguishes profile and label permutation spaces", {
  sizes <- rep(c(2L, 3L, 5L), each = 4L)
  unit_treatment <- rep(c("control", "control", "impact", "impact"), 3L)
  dat <- data.frame(
    site = factor(rep(seq_along(sizes), sizes)),
    treatment = factor(rep(unit_treatment, sizes))
  )
  set.seed(42)
  response <- stats::rnorm(length(sizes))[dat$site] +
    stats::rnorm(nrow(dat))
  design <- compile_permanova_design(
    ~ treatment + (1 | site), dat, dist(response)
  )
  fit <- permanova_unequal_random_intercept(
    design, permutations = 0
  )
  expect_equal(fit$permutation$attainable, 24 ^ 3)
  expect_equal(fit$permutation$attainable_label_assignments, 6 ^ 3)
  expect_equal(fit$test$attainable_label_assignments, 6 ^ 3)
})

test_that("audit separates valid nulls from sensitivity scenarios", {
  audit <- audit_unequal_random_intercept_calibration(
    batches = 1, simulations_per_batch = 5, permutations = 9,
    engine_checks = 0, seed = 20260719
  )
  validity <- unique(audit$summary[c("scenario", "exchangeability")])
  expect_true(all(validity$exchangeability[
    validity$scenario != "dispersion_violation"
  ]))
  expect_false(validity$exchangeability[
    validity$scenario == "dispersion_violation"
  ])
})

test_that("audit validates workloads and seed ranges", {
  expect_error(
    audit_unequal_random_intercept_calibration(
      batches = 1, simulations_per_batch = 2, engine_checks = 3
    ),
    "cannot exceed"
  )
  expect_error(
    audit_unequal_random_intercept_calibration(
      batches = 1, simulations_per_batch = 50000, engine_checks = 0
    ),
    "fewer than 50000"
  )
  expect_error(
    audit_unequal_random_intercept_calibration(
      batches = 3, simulations_per_batch = 2, engine_checks = 0,
      seed = .Machine$integer.max - 100
    ),
    "seed ranges"
  )
})
