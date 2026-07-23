.unbalanced_hybrid_example <- function(
    site_counts = c(R1 = 2L, R2 = 3L, R3 = 2L, R4 = 3L),
    missing_pairs = data.frame(treatment = "C", region = "R4"),
    treatment_specific_sizes = FALSE) {
  treatments <- c("A", "B", "C")
  rows <- list()
  index <- 0L
  for (region in names(site_counts)) {
    region_index <- match(region, names(site_counts))
    for (treatment in treatments) {
      omitted <- nrow(missing_pairs) && any(
        missing_pairs$treatment == treatment & missing_pairs$region == region
      )
      if (omitted) next
      for (site_index in seq_len(site_counts[[region]])) {
        size <- if (treatment_specific_sizes) {
          match(treatment, treatments)
        } else {
          1L + ((region_index + site_index) %% 2L)
        }
        index <- index + 1L
        rows[[index]] <- data.frame(
          treatment = treatment, region = region,
          site = paste(region, paste0("S", site_index), sep = "_"),
          replicate = seq_len(size), stringsAsFactors = FALSE
        )
      }
    }
  }
  data <- do.call(rbind, rows)
  data$treatment <- factor(data$treatment, levels = treatments)
  data$region <- factor(data$region, levels = names(site_counts))
  data$site <- factor(data$site)
  profile <- interaction(data$treatment, data$region, drop = TRUE)
  site_unit <- interaction(data$region, data$site, drop = TRUE)
  cell <- interaction(data$treatment, data$region, data$site, drop = TRUE)
  data$y <- c(A = -0.5, B = 0.3, C = 1.1)[data$treatment] +
    c(R1 = -0.4, R2 = 0.1, R3 = 0.5, R4 = -0.2)[data$region] +
    sin(as.integer(site_unit)) / 5 + cos(as.integer(profile)) / 4 +
    sin(as.integer(cell) + data$replicate) / 10
  data
}

.unbalanced_hybrid_formula <- function() {
  ~ treatment + (1 | region/site) +
    (1 | treatment:region) + (1 | treatment:region:site)
}

.profile_scalar_f <- function(response, structure) {
  y <- drop(structure$aggregated$membership %*% response)
  numerator <- sum(y * (structure$H_numerator %*% y)) /
    structure$numerator_df
  denominator <- sum(y * (structure$H_denominator %*% y)) /
    structure$denominator_df
  numerator / denominator
}

test_that("unbalanced hybrid equals the equal-site-weight profile GLH", {
  data <- .unbalanced_hybrid_example()
  design <- compile_permanova_design(
    .unbalanced_hybrid_formula(), data, dist(data$y)
  )
  fit <- permanova_unbalanced_hybrid(design, permutations = 0)
  structure <- distmer:::.validate_unbalanced_hybrid(design)

  expect_s3_class(fit, "permanova_unbalanced_hybrid")
  expect_equal(
    fit$test$statistic, .profile_scalar_f(data$y, structure),
    tolerance = 1e-10
  )
  expect_equal(fit$test$numerator$df, 2)
  expect_equal(fit$test$denominator$df, 5)
  expect_equal(fit$structure$observed_profiles, 11)
  expect_match(fit$estimand, "equal.*site")
  expect_match(fit$test$warnings[[1]], "no unbalanced parametric EMS")
})

test_that("unbalanced hybrid exact permutations equal direct profile responses", {
  data <- .unbalanced_hybrid_example()
  design <- compile_permanova_design(
    .unbalanced_hybrid_formula(), data, dist(data$y)
  )
  fit <- permanova_unbalanced_hybrid(
    design, permutations = "all", keep_permutations = TRUE,
    max_exact = 1000
  )
  structure <- distmer:::.validate_unbalanced_hybrid(design)
  plan <- distmer:::.stratified_permutation_plan(
    structure$exchangeability_strata, "all", 1000
  )
  profile_y <- drop(structure$aggregated$membership %*% data$y)
  direct <- apply(plan$orders, 1L, function(order) {
    y <- profile_y[order]
    (sum(y * (structure$H_numerator %*% y)) / structure$numerator_df) /
      (sum(y * (structure$H_denominator %*% y)) /
         structure$denominator_df)
  })

  expect_true(fit$test$exact_permutation)
  expect_equal(fit$test$attainable_permutations, factorial(3)^3 * factorial(2))
  expect_equal(fit$test$permutation_statistics, direct, tolerance = 1e-10)
})

test_that("unbalanced hybrid routes through mixed contrasts and max-F", {
  data <- .unbalanced_hybrid_example()
  design <- compile_permanova_design(
    .unbalanced_hybrid_formula(), data, dist(data$y)
  )
  contrast <- c(A = 1, B = -1, C = 0)
  ordinary <- permanova_mixed_contrasts(
    design, contrast, permutations = 0
  )
  maximum <- permanova_mixed_maxT(
    design, contrast, permutations = "all", max_exact = 1000
  )
  structure <- distmer:::.mixed_contrast_structure(design)
  units <- distmer:::.prepare_mixed_contrast_units(structure)
  y <- drop(structure$aggregated$membership %*% data$y)
  coefficient_row <- as.numeric(contrast %*% units$mean_rows)
  projector <- distmer:::.linear_contrast_projector(
    units$X, coefficient_row, units$H_full, 1e-10
  )
  direct <- sum(y * (projector$H %*% y)) /
    (sum(y * (units$H_denominator %*% y)) / units$denominator_df)

  expect_match(ordinary$design_type, "unbalanced crossed-nested")
  expect_equal(ordinary$table$pseudo_f, direct, tolerance = 1e-10)
  expect_equal(maximum$table$pseudo_f, direct, tolerance = 1e-10)
  expect_equal(maximum$table$p_perm, maximum$table$p_maxT)
  expect_match(ordinary$estimand, "additive region adjustment")
})

test_that("unbalanced hybrid rejects unsafe incidence and topology", {
  balanced <- .unbalanced_hybrid_example(
    site_counts = c(R1 = 2L, R2 = 2L, R3 = 2L, R4 = 2L),
    missing_pairs = data.frame(treatment = character(), region = character())
  )
  balanced <- balanced[balanced$replicate == 1L, , drop = FALSE]
  balanced_design <- compile_permanova_design(
    .unbalanced_hybrid_formula(), balanced, dist(balanced$y)
  )
  expect_error(
    permanova_unbalanced_hybrid(balanced_design, permutations = 0),
    "balanced; use `permanova_hybrid"
  )

  partial <- .unbalanced_hybrid_example()
  partial <- partial[!(
    partial$treatment == "B" & partial$region == "R1" &
      partial$site == "R1_S1"
  ), ]
  partial_design <- compile_permanova_design(
    .unbalanced_hybrid_formula(), partial, dist(partial$y)
  )
  expect_error(
    permanova_unbalanced_hybrid(partial_design, permutations = 0),
    "complete site set"
  )

  confounded <- .unbalanced_hybrid_example(treatment_specific_sizes = TRUE)
  confounded_design <- compile_permanova_design(
    .unbalanced_hybrid_formula(), confounded, dist(confounded$y)
  )
  expect_error(
    permanova_unbalanced_hybrid(confounded_design, permutations = 0),
    "topology strata do not span"
  )

  disconnected_pairs <- expand.grid(
    treatment = c("C", "A", "B"), region = c("R1", "R2", "R3", "R4"),
    stringsAsFactors = FALSE
  )
  disconnected_pairs <- disconnected_pairs[
    (disconnected_pairs$treatment == "C" &
       disconnected_pairs$region %in% c("R1", "R2")) |
      (disconnected_pairs$treatment %in% c("A", "B") &
         disconnected_pairs$region %in% c("R3", "R4")), , drop = FALSE
  ]
  disconnected <- .unbalanced_hybrid_example(
    missing_pairs = disconnected_pairs
  )
  disconnected_design <- compile_permanova_design(
    .unbalanced_hybrid_formula(), disconnected, dist(disconnected$y)
  )
  expect_error(
    permanova_unbalanced_hybrid(disconnected_design, permutations = 0),
    "incidence is disconnected"
  )
})

test_that("unbalanced hybrid rejects zero denominators and retains semimetrics", {
  data <- .unbalanced_hybrid_example()
  additive <- c(A = -1, B = 0, C = 1)[data$treatment] +
    c(R1 = -.3, R2 = .2, R3 = .5, R4 = -.1)[data$region]
  zero_design <- compile_permanova_design(
    .unbalanced_hybrid_formula(), data, dist(additive)
  )
  expect_error(
    permanova_unbalanced_hybrid(zero_design, permutations = 0),
    "residual mean square is numerically zero"
  )

  distance <- as.matrix(dist(cbind(data$y, seq_len(nrow(data)) / 10)))
  distance[1, nrow(distance)] <- distance[nrow(distance), 1] <-
    5 * max(distance)
  semimetric <- compile_permanova_design(
    .unbalanced_hybrid_formula(), data, distance
  )
  fit <- permanova_unbalanced_hybrid(
    semimetric, permutations = 19, seed = 26
  )
  expect_false(semimetric$spectral$euclidean)
  expect_true(is.finite(fit$test$p_value))
  expect_match(paste(fit$test$warnings, collapse = " "), "non-Euclidean")
})
