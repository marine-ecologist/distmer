.direct_pairwise_statistics <- function(design, contrasts, orders) {
  structure <- distmer:::.validate_contrast_design(design)
  projectors <- lapply(seq_len(nrow(contrasts)), function(index) {
    distmer:::.contrast_projector(
      contrasts[index, ], structure$group, structure$counts,
      structure$H_full, 1e-10
    )$H
  })
  t(vapply(seq_len(nrow(orders)), function(row) {
    G <- design$G[orders[row, ], orders[row, ], drop = FALSE]
    denominator <- distmer:::trace_hgh(
      structure$H_denominator, G
    ) / structure$denominator_df
    vapply(projectors, function(H) {
      distmer:::trace_hgh(H, G) / denominator
    }, numeric(1))
  }, numeric(nrow(contrasts))))
}

test_that("closed pairwise exact intersections match direct reconstruction", {
  data <- data.frame(
    group = factor(rep(c("A", "B", "C"), each = 2)),
    y = c(-1, 0, .2, 1.1, 1.7, 3)
  )
  design <- compile_permanova_design(~ group, data, dist(data$y))
  fit <- permanova_pairwise_strong(
    design, permutations = "all", keep_permutations = TRUE
  )

  expect_s3_class(fit, "permanova_pairwise_strong")
  expect_equal(nrow(fit$intersections), 7L)
  expect_equal(fit$intersections$attainable[1:3], rep(6, 3))
  expect_equal(fit$intersections$attainable[4:7], rep(90, 4))
  for (detail in fit$intersection_details) {
    direct <- .direct_pairwise_statistics(
      design, fit$contrast_matrix, detail$orders
    )
    direct_maximum <- apply(
      direct[, detail$subset, drop = FALSE], 1L, max
    )
    expect_equal(detail$statistics, direct_maximum, tolerance = 1e-12)
    expect_equal(
      detail$p_value,
      mean(direct_maximum >= detail$statistic - 1e-10)
    )
  }
  for (index in seq_len(nrow(fit$table))) {
    containing <- vapply(fit$intersection_details, function(detail) {
      index %in% detail$subset
    }, logical(1))
    expect_equal(
      fit$table$p_closed[[index]],
      max(fit$intersections$p_perm[containing])
    )
  }
  expect_true(all(fit$table$p_closed >= fit$table$p_perm))
})

test_that("partial-null permutations keep non-null components fixed", {
  data <- data.frame(
    group = factor(rep(c("A", "B", "C", "D"), each = 2)),
    y = c(-1, 1, -1, 1, 9, 11, 19, 21)
  )
  design <- compile_permanova_design(~ group, data, dist(data$y))
  fit <- permanova_pairwise_strong(
    design, permutations = "all", keep_permutations = TRUE,
    max_exact = 10000
  )
  ab <- fit$intersection_details[[which(vapply(
    fit$intersection_details,
    function(detail) identical(detail$hypotheses, "A - B"), logical(1)
  ))]]
  expect_true(all(ab$orders[, 5:8, drop = FALSE] ==
                    matrix(5:8, nrow(ab$orders), 4, byrow = TRUE)))
  expect_equal(ab$attainable, 6)

  disconnected <- fit$intersection_details[[which(vapply(
    fit$intersection_details, function(detail) {
      setequal(detail$hypotheses, c("A - B", "C - D"))
    }, logical(1)))]]
  expect_equal(disconnected$attainable, 36)
  expect_true(all(disconnected$orders[, 1:4] %in% 1:4))
  expect_true(all(disconnected$orders[, 5:8] %in% 5:8))
  expect_match(fit$exchangeability_proof, "non-null groups", fixed = TRUE,
               all = FALSE)
  expect_match(fit$error_control, "strong FWER")
})

test_that("Monte Carlo closed tests are reproducible", {
  data <- data.frame(
    group = factor(rep(c("A", "B", "C"), each = 4)),
    y = c(-2, -.4, .2, 1, -1, -.2, .4, 1.4, 1, 2, 3, 4)
  )
  design <- compile_permanova_design(~ group, data, dist(data$y))
  first <- permanova_pairwise_strong(
    design, permutations = 29, seed = 18, keep_permutations = TRUE
  )
  second <- permanova_pairwise_strong(
    design, permutations = 29, seed = 18, keep_permutations = TRUE
  )
  expect_equal(first$table, second$table)
  expect_equal(
    lapply(first$intersection_details, `[[`, "orders"),
    lapply(second$intersection_details, `[[`, "orders")
  )
  expect_false(first$permutation$exact)
  expect_equal(first$permutation$total_used, 7L * 29L)
})

test_that("strong pairwise inference enforces its proof boundary", {
  data <- data.frame(
    group = factor(rep(c("A", "B", "C"), each = 2)),
    block = factor(rep(1:2, 3)),
    y = c(-1, 0, .2, 1.1, 1.7, 3)
  )
  mixed <- compile_permanova_design(
    ~ group + (1 | block), data, dist(data$y)
  )
  expect_error(permanova_pairwise_strong(mixed), "fixed one-factor")

  multifactor <- compile_permanova_design(
    ~ group + block, data, dist(data$y)
  )
  expect_error(permanova_pairwise_strong(multifactor), "exactly one")

  many <- data.frame(group = factor(rep(LETTERS[1:6], each = 2)), y = 1:12)
  many_design <- compile_permanova_design(~ group, many, dist(many$y))
  expect_error(
    permanova_pairwise_strong(many_design, permutations = 0),
    "exceeding `max_intersections"
  )
  expect_error(
    permanova_pairwise_strong(
      compile_permanova_design(~ group, data, dist(data$y)),
      permutations = "all", max_exact = 10
    ),
    "exceeds `max_exact"
  )
  expect_error(
    permanova_pairwise_strong(
      compile_permanova_design(~ group, data, dist(data$y)),
      permutations = 1, max_exact = 0
    ),
    "positive integer"
  )

  constant <- data.frame(
    group = factor(rep(c("A", "B"), each = 2)), y = rep(1, 4)
  )
  expect_error(
    permanova_pairwise_strong(
      compile_permanova_design(~ group, constant, dist(constant$y)),
      permutations = 0
    ),
    "numerically zero"
  )
})

test_that("strong closed tests retain non-Euclidean geometry", {
  data <- data.frame(
    group = factor(rep(c("A", "B", "C"), each = 2)),
    y = c(-1, 0, .2, 1.1, 1.7, 3)
  )
  semimetric <- as.matrix(dist(data$y))
  semimetric[1, 6] <- semimetric[6, 1] <- 10
  design <- compile_permanova_design(~ group, data, semimetric)
  fit <- permanova_pairwise_strong(design, permutations = 9, seed = 2)
  expect_false(design$spectral$euclidean)
  expect_true(all(is.finite(fit$table$p_closed)))
  expect_match(fit$tests[[1]]$warnings[[3]], "non-Euclidean")
})
