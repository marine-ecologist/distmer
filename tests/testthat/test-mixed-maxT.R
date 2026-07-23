.mixed_maxT_block_example <- function() {
  data <- expand.grid(
    replicate = factor(1:2),
    treatment = factor(c("A", "B", "C"), levels = c("A", "B", "C")),
    block = factor(1:2)
  )
  data$y <- c(A = -0.4, B = 0.3, C = 1.1)[data$treatment] +
    c(`1` = -0.2, `2` = 0.5)[data$block] +
    c(-0.2, 0.1, 0.1, 0.3, -0.1, -0.2)[
      interaction(data$treatment, data$block)
    ] +
    rep(c(-0.15, 0.15), length.out = nrow(data))
  compile_permanova_design(
    ~ treatment + (1 | block) + (1 | treatment:block),
    data, dist(data$y)
  )
}

test_that("mixed max-F equals direct shared unit permutations", {
  design <- .mixed_maxT_block_example()
  contrasts <- rbind(
    `A - B` = c(A = 1, B = -1, C = 0),
    `B - C` = c(A = 0, B = 1, C = -1)
  )
  fit <- permanova_mixed_maxT(
    design, contrasts, permutations = "all",
    keep_permutations = TRUE, max_exact = 1000
  )
  structure <- distmer:::.mixed_contrast_structure(design)
  units <- distmer:::.prepare_mixed_contrast_units(structure)
  plan <- distmer:::.mixed_contrast_plan(
    structure, "all", 1000
  )
  projectors <- lapply(seq_len(nrow(contrasts)), function(index) {
    row <- as.numeric(contrasts[index, ] %*% units$mean_rows)
    distmer:::.linear_contrast_projector(
      units$X, row, units$H_full, 1e-10
    )$H
  })
  direct <- t(apply(plan$orders, 1L, function(order) {
    G <- structure$aggregated$G[order, order, drop = FALSE]
    denominator <- distmer:::trace_hgh(
      units$H_denominator, G
    ) / units$denominator_df
    vapply(projectors, function(H) {
      distmer:::trace_hgh(H, G) / denominator
    }, numeric(1))
  }))

  expect_equal(fit$permutation$used, 36)
  expect_true(fit$permutation$exact)
  expect_equal(unname(fit$permutation$statistic_matrix), direct,
               tolerance = 1e-10)
  expect_equal(fit$permutation$maxima, apply(direct, 1L, max),
               tolerance = 1e-10)
  expect_equal(fit$table$p_maxT, vapply(
    fit$table$pseudo_f,
    function(observed) mean(apply(direct, 1L, max) >= observed - 1e-10),
    numeric(1)
  ))
})

test_that("one mixed contrast has identical raw and maximum p-values", {
  design <- .mixed_maxT_block_example()
  contrast <- c(A = 1, B = -1, C = 0)
  fit <- permanova_mixed_maxT(
    design, contrast, permutations = "all", max_exact = 1000
  )
  ordinary <- permanova_mixed_contrasts(
    design, contrast, permutations = 0
  )

  expect_equal(fit$table$p_perm, fit$table$p_maxT)
  expect_equal(fit$table$pseudo_f, ordinary$table$pseudo_f,
               tolerance = 1e-12)
})

test_that("mixed pairwise max-F retains exact-size unit strata", {
  sizes <- rep(2:4, each = 3)
  treatment_by_site <- factor(
    rep(c("A", "B", "C"), 3), levels = c("A", "B", "C")
  )
  data <- data.frame(
    site = factor(rep(seq_along(sizes), sizes)),
    treatment = factor(rep(treatment_by_site, sizes),
                       levels = levels(treatment_by_site))
  )
  data$y <- rep(seq_along(sizes) / 5, sizes) +
    c(A = -0.3, B = 0.2, C = 0.7)[data$treatment] +
    sin(seq_len(nrow(data))) / 10
  design <- compile_permanova_design(
    ~ treatment + (1 | site), data, dist(data$y)
  )
  fit <- permanova_mixed_maxT_pairwise(
    design, permutations = 19, seed = 4, keep_permutations = TRUE
  )

  expect_equal(nrow(fit$table), 3)
  expect_true(all(fit$table$p_maxT >= fit$table$p_perm))
  expect_match(fit$permutation$algorithm, "exact random-intercept unit")
  expect_equal(ncol(fit$permutation$statistic_matrix), 3)
  expect_match(fit$error_control, "weak familywise")
  expect_match(fit$error_control, "strong control is not claimed")
})

test_that("mixed max-F rejects unsupported structures and retains semimetric geometry", {
  design <- .mixed_maxT_block_example()
  distance <- as.matrix(design$distance)
  distance[1, nrow(distance)] <- distance[nrow(distance), 1] <-
    5 * max(distance)
  semimetric <- compile_permanova_design(
    ~ treatment + (1 | block) + (1 | treatment:block),
    design$data, distance
  )
  fit <- permanova_mixed_maxT(
    semimetric, c(A = 1, B = -1, C = 0),
    permutations = 9, seed = 5
  )
  expect_false(semimetric$spectral$euclidean)
  expect_true(is.finite(fit$table$p_maxT))
  expect_match(paste(fit$tests[[1]]$warnings, collapse = " "),
               "non-Euclidean")

  crossed_data <- expand.grid(
    treatment = factor(c("A", "B", "C")),
    block = factor(1:2), observer = factor(1:2)
  )
  crossed_data$y <- seq_len(nrow(crossed_data))
  unsupported <- compile_permanova_design(
    ~ treatment + (1 | block) + (1 | observer),
    crossed_data, dist(crossed_data$y)
  )
  expect_error(
    permanova_mixed_maxT(
      unsupported, c(A = 1, B = -1, C = 0), permutations = 0
    ),
    "support one random intercept"
  )
})
