.manual_maxT <- function(design, contrasts) {
  structure <- distmer:::.validate_contrast_design(design)
  plan <- distmer:::.permutation_orders(design$n, "all", 100000)
  projectors <- lapply(seq_len(nrow(contrasts)), function(index) {
    distmer:::.contrast_projector(
      contrasts[index, ], structure$group, structure$counts,
      structure$H_full, 1e-10
    )$H
  })
  calculate <- function(G) {
    denominator <- distmer:::trace_hgh(
      structure$H_denominator, G
    ) / structure$denominator_df
    vapply(projectors, function(H) {
      distmer:::trace_hgh(H, G) / denominator
    }, numeric(1))
  }
  observed <- calculate(design$G)
  statistics <- t(apply(plan$orders, 1, function(order) {
    calculate(design$G[order, order, drop = FALSE])
  }))
  maxima <- apply(statistics, 1, max)
  list(
    observed = observed, statistics = statistics, maxima = maxima,
    raw = vapply(seq_along(observed), function(index) {
      mean(statistics[, index] >= observed[[index]] - 1e-10)
    }, numeric(1)),
    adjusted = vapply(observed, function(value) {
      mean(maxima >= value - 1e-10)
    }, numeric(1))
  )
}

test_that("exact max-F equals independent brute-force enumeration", {
  data <- data.frame(
    group = factor(rep(c("A", "B", "C"), each = 2)),
    y = c(-1, 0, .2, 1.1, 1.7, 3)
  )
  design <- compile_permanova_design(~ group, data, dist(data$y))
  contrasts <- rbind(
    `A - B` = c(A = 1, B = -1, C = 0),
    `AB - C` = c(A = .5, B = .5, C = -1)
  )
  fit <- permanova_maxT(
    design, contrasts, permutations = "all", keep_permutations = TRUE
  )
  reference <- .manual_maxT(design, contrasts)

  expect_equal(fit$table$pseudo_f, reference$observed, tolerance = 1e-12)
  expect_equal(unname(fit$permutation$statistic_matrix), reference$statistics,
               tolerance = 1e-12)
  expect_equal(fit$permutation$maxima, reference$maxima, tolerance = 1e-12)
  expect_equal(fit$table$p_perm, reference$raw)
  expect_equal(fit$table$p_maxT, reference$adjusted)
  expect_true(all(fit$table$p_maxT >= fit$table$p_perm))
  expect_equal(fit$permutation$used, factorial(6))
})

test_that("one contrast has identical raw and max-F p-values", {
  data <- data.frame(
    group = factor(rep(c("A", "B", "C"), each = 2)),
    y = c(-1, 0, .2, 1.1, 1.7, 3)
  )
  design <- compile_permanova_design(~ group, data, dist(data$y))
  fit <- permanova_maxT(
    design, c(A = 1, B = -1, C = 0), permutations = "all"
  )
  expect_equal(fit$table$p_perm, fit$table$p_maxT)
  expect_match(fit$error_control, "weak familywise")
  expect_match(fit$tests[[1]]$warnings[[1]], "strong subset-null control is not claimed")
})

test_that("pairwise max-F uses one shared family distribution", {
  data <- data.frame(
    group = factor(rep(c("A", "B", "C"), c(3, 2, 3))),
    y = c(-1, 0, .1, .4, 1, 1.3, 2, 3)
  )
  design <- compile_permanova_design(~ group, data, dist(data$y))
  first <- permanova_maxT_pairwise(
    design, permutations = 29, seed = 8, keep_permutations = TRUE
  )
  second <- permanova_maxT_pairwise(
    design, permutations = 29, seed = 8, keep_permutations = TRUE
  )

  expect_equal(nrow(first$table), 3)
  expect_equal(first$permutation$statistic_matrix,
               second$permutation$statistic_matrix)
  expect_equal(first$table$p_maxT, second$table$p_maxT)
  expect_equal(first$permutation$maxima,
               apply(first$permutation$statistic_matrix, 1, max))
})

test_that("max-F rejects unsupported models and retains semimetric geometry", {
  data <- data.frame(
    group = factor(rep(c("A", "B", "C"), each = 2)),
    block = factor(rep(1:2, 3)), y = c(-1, 0, .2, 1.1, 1.7, 3)
  )
  mixed <- compile_permanova_design(
    ~ group + (1 | block), data, dist(data$y)
  )
  expect_error(permanova_maxT_pairwise(mixed), "fixed one-factor")

  multifactor <- compile_permanova_design(
    ~ group + block, data, dist(data$y)
  )
  expect_error(permanova_maxT_pairwise(multifactor), "exactly one")

  semimetric <- as.matrix(dist(data$y))
  semimetric[1, 6] <- semimetric[6, 1] <- 10
  design <- compile_permanova_design(~ group, data, semimetric)
  fit <- permanova_maxT_pairwise(design, permutations = 9, seed = 2)
  expect_false(design$spectral$euclidean)
  expect_true(all(is.finite(fit$table$p_maxT)))
  expect_match(fit$tests[[1]]$warnings[[2]], "non-Euclidean")
})
