test_that("balanced simulation is reproducible and correctly shaped", {
  first <- simulate_balanced_mixed(
    groups_per_treatment = 3, replicates = 2, response_variables = 4, seed = 10
  )
  second <- simulate_balanced_mixed(
    groups_per_treatment = 3, replicates = 2, response_variables = 4, seed = 10
  )
  expect_equal(first$response, second$response)
  expect_equal(dim(first$response), c(12, 4))
  expect_equal(as.integer(table(first$data$site)), rep(2L, 6))
})

test_that("deterministic validation suite passes", {
  report <- run_permanova_validation(seed = 4)
  expect_s3_class(report, "permanova_validation")
  expect_true(all(report$passed), info = paste(report$detail[!report$passed], collapse = "; "))
})

test_that("simulation framework reports all operating characteristics", {
  report <- simulate_mixed_validation(
    simulations = 2, permutations = 3, groups_per_treatment = 3,
    replicates = 2, seed = 4, keep_results = TRUE
  )
  expect_s3_class(report, "permanova_simulation_validation")
  expect_equal(report$summary$scenario,
               c("null", "location_shift", "dispersion_heterogeneity"))
  expect_equal(nrow(report$results), 6)
  expect_true(all(report$summary$failures == 0))
  expect_lt(max(report$summary$max_abs_denominator_f_difference), 1e-10)
})

test_that("benchmark harness returns reproducible schema", {
  timings <- benchmark_permanova(
    group_counts = 4, replicates = 2, permutations = 3, seed = 2
  )
  expect_named(timings, c(
    "observations", "whole_plots", "replicates", "permutations", "elapsed_seconds"
  ))
  expect_equal(timings$observations, 8)
  expect_gte(timings$elapsed_seconds, 0)
})

test_that("mixed contrast validation matches the scalar GLH", {
  validation <- simulate_mixed_contrast_validation(
    simulations = 3, permutations = 7, blocks = 4,
    replicates = 2, seed = 20260719
  )
  expect_s3_class(validation, "permanova_mixed_contrast_validation")
  expect_lt(validation$max_engine_classical_f_difference, 1e-10)
  expect_equal(nrow(validation$summary), 4)
  expect_true(all(is.finite(validation$summary$confidence_low)))
})

test_that("max-F validation reports complete-null family error", {
  validation <- simulate_maxT_validation(
    simulations = 3, permutations = 7, seed = 20260719
  )
  expect_s3_class(validation, "permanova_maxT_validation")
  expect_equal(nrow(validation$summary), 6)
  expect_equal(unique(validation$summary$method),
               c("unadjusted", "holm", "maxT"))
  expect_match(validation$interpretation, "weak-FWER")
})

test_that("unbalanced mixed contrast validation matches direct references", {
  validation <- simulate_unbalanced_mixed_contrast_validation(
    simulations = 2, permutations = 7, seed = 20260720
  )
  expect_s3_class(
    validation, "permanova_unbalanced_mixed_contrast_validation"
  )
  expect_lt(validation$max_engine_classical_f_difference, 1e-10)
  expect_equal(validation$max_engine_direct_p_difference, 0)
  expect_lt(validation$max_engine_direct_permutation_difference, 1e-10)
})

test_that("declared crossed interaction agrees with the classical test", {
  validation <- simulate_crossed_interaction_validation(
    simulations = 100, blocks = 8, replicates = 3,
    interaction_sd = 2, seed = 20260718
  )
  expect_s3_class(validation, "permanova_crossed_validation")
  expect_lt(validation$max_declared_classical_f_difference, 1e-10)
  expect_lt(validation$max_declared_classical_p_difference, 1e-10)
  expect_equal(
    validation$summary$rejection_rate[[1L]],
    validation$summary$rejection_rate[[2L]]
  )
  expect_gt(
    validation$summary$rejection_rate[[3L]],
    validation$summary$rejection_rate[[1L]] + 0.2
  )
  expect_true(all(is.finite(validation$summary$confidence_low)))
  expect_true(all(is.finite(validation$summary$confidence_high)))
})

test_that("nested validation agrees with region-mean ANOVA", {
  validation <- simulate_nested_validation(
    simulations = 100, regions_per_treatment = 4,
    sites_per_region = 2, replicates = 3, seed = 20260718
  )
  expect_s3_class(validation, "permanova_nested_validation")
  expect_lt(validation$max_nested_classical_f_difference, 1e-10)
  expect_lt(validation$max_nested_classical_p_difference, 1e-10)
  expect_equal(validation$summary$rejection_rate[[1L]],
               validation$summary$rejection_rate[[2L]])
  expect_true(all(is.finite(validation$summary$confidence_low)))
  expect_true(all(is.finite(validation$summary$confidence_high)))
})

test_that("level-specific nested validation agrees at every hierarchy level", {
  validation <- simulate_nested_level_validation(
    simulations = 100, regions_per_treatment = 4,
    replicates = 3, seed = 20260718
  )
  expect_s3_class(validation, "permanova_nested_level_validation")
  expect_lt(max(validation$max_f_difference), 1e-10)
  expect_lt(max(validation$max_p_difference), 1e-10)
  for (level in c("region", "site", "observation")) {
    rates <- validation$summary$rejection_rate[
      validation$summary$level == level
    ]
    expect_equal(rates[[1L]], rates[[2L]])
  }
  expect_true(all(is.finite(validation$summary$confidence_low)))
  expect_true(all(is.finite(validation$summary$confidence_high)))
})

test_that("crossed permutation validation matches the direct reference", {
  validation <- simulate_crossed_permutation_validation(
    simulations = 20, permutations = 19, blocks = 8,
    replicates = 2, seed = 20260718
  )
  expect_s3_class(
    validation, "permanova_crossed_permutation_validation"
  )
  expect_lt(validation$max_engine_direct_f_difference, 1e-10)
  expect_equal(validation$max_engine_direct_p_difference, 0)
  expect_lt(
    validation$max_engine_direct_permutation_statistic_difference,
    1e-10
  )
  expect_lt(validation$max_engine_classical_f_difference, 1e-10)
  for (scenario in unique(validation$summary$scenario)) {
    values <- validation$summary[
      validation$summary$scenario == scenario, , drop = FALSE
    ]
    expect_equal(values$rejection_rate[[1L]],
                 values$rejection_rate[[2L]])
  }
  expect_true(all(is.finite(validation$summary$confidence_low)))
  expect_true(all(is.finite(validation$summary$confidence_high)))
})

test_that("dispersion validation matches scalar references across structures", {
  validation <- simulate_dispersion_validation(
    simulations = 2, permutations = 3, seed = 20260718
  )
  expect_s3_class(validation, "permanova_dispersion_validation")
  expect_equal(nrow(validation$summary), 12)
  expect_lt(validation$max_engine_classical_f_difference, 1e-10)
  expect_true(all(is.finite(validation$summary$confidence_low)))
  expect_true(all(is.finite(validation$summary$confidence_high)))
})

test_that("hybrid validation matches direct treatment-profile inference", {
  validation <- simulate_hybrid_validation(
    simulations = 5, permutations = 7, regions = 4,
    sites = 2, replicates = 2, seed = 20260718
  )
  expect_s3_class(validation, "permanova_hybrid_validation")
  expect_lt(validation$max_engine_direct_f_difference, 1e-10)
  expect_equal(validation$max_engine_direct_p_difference, 0)
  expect_lt(
    validation$max_engine_direct_permutation_statistic_difference,
    1e-10
  )
  expect_lt(validation$max_engine_classical_f_difference, 1e-10)
  for (scenario in unique(validation$summary$scenario)) {
    values <- validation$summary[
      validation$summary$scenario == scenario, , drop = FALSE
    ]
    expect_equal(values$rejection_rate[[1L]],
                 values$rejection_rate[[2L]])
  }
})

test_that("unbalanced validation matches direct Freedman-Lane inference", {
  validation <- simulate_unbalanced_validation(
    simulations = 5, permutations = 7, seed = 20260718
  )
  expect_s3_class(validation, "permanova_unbalanced_validation")
  expect_lt(validation$max_engine_direct_f_difference, 1e-10)
  expect_equal(validation$max_engine_direct_p_difference, 0)
  expect_lt(
    validation$max_engine_direct_permutation_statistic_difference,
    1e-10
  )
  expect_lt(validation$max_engine_classical_f_difference, 1e-10)
  for (scenario in unique(validation$summary$scenario)) {
    values <- validation$summary[
      validation$summary$scenario == scenario, , drop = FALSE
    ]
    expect_equal(values$rejection_rate[[1L]],
                 values$rejection_rate[[2L]])
  }
})

test_that("unbalanced mixed validation matches unit-level references", {
  validation <- simulate_unbalanced_mixed_validation(
    simulations = 3, permutations = 7, seed = 20260719
  )
  expect_s3_class(validation, "permanova_unbalanced_mixed_validation")
  expect_lt(validation$max_engine_direct_f_difference, 1e-10)
  expect_equal(validation$max_engine_direct_p_difference, 0)
  expect_lt(
    validation$max_engine_direct_permutation_statistic_difference,
    1e-10
  )
  for (structure in c("nested", "block")) {
    for (scenario in c("exchangeable_null", "location_shift")) {
      rates <- validation$summary$rejection_rate[
        validation$summary$structure == structure &
          validation$summary$scenario == scenario
      ]
      expect_equal(rates[[1L]], rates[[2L]])
    }
  }
})

test_that("planned contrast validation matches constrained scalar inference", {
  validation <- simulate_contrast_validation(
    simulations = 3, permutations = 7, seed = 20260719
  )
  expect_s3_class(validation, "permanova_contrast_validation")
  expect_lt(validation$max_engine_direct_f_difference, 1e-10)
  expect_equal(validation$max_engine_direct_p_difference, 0)
  expect_lt(
    validation$max_engine_direct_permutation_statistic_difference,
    1e-10
  )
  expect_lt(validation$max_engine_classical_f_difference, 1e-10)
  for (scenario in c("exchangeable_null", "location_shift")) {
    rates <- validation$summary$rejection_rate[
      validation$summary$scenario == scenario &
        validation$summary$method %in% c("engine", "direct")
    ]
    expect_equal(rates[[1L]], rates[[2L]])
  }
})

test_that("marginal contrast validation matches scalar GLH inference", {
  validation <- simulate_marginal_contrast_validation(
    simulations = 3, permutations = 7, seed = 20260719
  )
  expect_s3_class(validation, "permanova_marginal_contrast_validation")
  expect_lt(validation$max_engine_direct_f_difference, 1e-10)
  expect_equal(validation$max_engine_direct_p_difference, 0)
  expect_lt(
    validation$max_engine_direct_permutation_statistic_difference,
    1e-10
  )
  expect_lt(validation$max_engine_classical_f_difference, 1e-10)
})

test_that("adjusted contrast validation matches scalar and direct inference", {
  validation <- simulate_adjusted_contrast_validation(
    simulations = 3, permutations = 7, seed = 20260719
  )
  expect_s3_class(validation, "permanova_adjusted_contrast_validation")
  expect_lt(validation$max_engine_classical_f_difference, 1e-10)
  expect_equal(validation$max_engine_direct_p_difference, 0)
  expect_lt(validation$max_engine_direct_permutation_difference, 1e-10)
  for (scenario in c("null", "location_shift")) {
    rates <- validation$summary$rejection_rate[
      validation$summary$scenario == scenario &
        validation$summary$method %in% c("engine", "direct")
    ]
    expect_equal(rates[[1L]], rates[[2L]])
  }
})

test_that("mixed max-F validation matches scalar unit permutations", {
  validation <- simulate_mixed_maxT_validation(
    simulations = 3, permutations = 7, seed = 20260719
  )
  expect_s3_class(validation, "permanova_mixed_maxT_validation")
  expect_lt(validation$max_engine_scalar_f_difference, 1e-10)
  expect_lt(validation$max_engine_scalar_statistic_difference, 1e-10)
  expect_equal(validation$max_engine_scalar_p_difference, 0)
})

test_that("lower unbalanced nested validation matches scalar references", {
  validation <- simulate_unbalanced_nested_lower_validation(
    simulations = 2, permutations = 7, seed = 20260719
  )
  expect_s3_class(
    validation, "permanova_unbalanced_nested_lower_validation"
  )
  expect_lt(validation$max_engine_classical_f_difference, 1e-10)
  expect_equal(validation$max_engine_direct_p_difference, 0)
  expect_lt(validation$max_engine_direct_permutation_difference, 1e-10)
})

test_that("unbalanced hybrid validation matches profile references", {
  validation <- simulate_unbalanced_hybrid_validation(
    simulations = 2, permutations = 7, seed = 20260724
  )
  expect_s3_class(
    validation, "permanova_unbalanced_hybrid_validation"
  )
  expect_lt(validation$max_engine_classical_f_difference, 1e-10)
  expect_equal(validation$max_engine_direct_p_difference, 0)
  expect_lt(validation$max_engine_direct_permutation_difference, 1e-10)
})

test_that("random-slope validation matches scalar and sign-flip references", {
  validation <- simulate_random_slope_validation(
    simulations = 2, permutations = 7, seed = 20260725
  )
  expect_s3_class(validation, "permanova_random_slope_validation")
  expect_lt(validation$max_engine_scalar_f_difference, 1e-10)
  expect_equal(validation$max_engine_direct_p_difference, 0)
  expect_lt(validation$max_engine_direct_permutation_difference, 1e-10)
  expect_lt(validation$max_engine_scalar_covariance_difference, 1e-10)
  expect_match(validation$interpretation, "central")
})

test_that("strong pairwise validation preserves partial-null components", {
  validation <- simulate_strong_pairwise_validation(
    simulations = 2, permutations = 7, seed = 20260726
  )
  expect_s3_class(validation, "permanova_strong_pairwise_validation")
  expect_lt(validation$max_direct_statistic_difference, 1e-10)
  expect_equal(validation$max_direct_p_difference, 0)
  expect_equal(validation$component_preservation_violations, 0)
  expect_match(validation$proof_boundary, "complete distributions")
})
