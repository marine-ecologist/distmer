test_that("quasi denominator matches a target EMS", {
  candidates <- rbind(low = c(residual = 1, random = 0),
                      high = c(residual = 1, random = 4))
  result <- quasi_denominator(
    c(residual = 1, random = 2), candidates,
    mean_squares = c(3, 7), dfs = c(10, 12)
  )

  expect_equal(result$coefficients, c(low = 0.5, high = 0.5), tolerance = 1e-10)
  expect_equal(result$ems_achieved, c(residual = 1, random = 2), tolerance = 1e-10)
  expect_equal(result$value, 5, tolerance = 1e-10)
  expect_false(result$exact_ems_match)
  expect_false(result$exact_single_mean_square)
  expect_gt(result$df, 0)
})

test_that("matching orthogonal denominator strata are pooled", {
  candidates <- rbind(a = c(residual = 1, random = 0),
                      b = c(residual = 1, random = 0))
  result <- quasi_denominator(
    c(residual = 1, random = 0), candidates,
    mean_squares = c(2, 4), dfs = c(3, 9)
  )
  expect_true(result$exact_ems_match)
  expect_equal(result$coefficients, c(a = 0.25, b = 0.75))
  expect_equal(result$value, 3.5)
  expect_equal(result$df, 12)
})

.crossed_example <- function() {
  dat <- expand.grid(
    replicate = seq_len(2), treatment = factor(c("control", "impact")),
    block = factor(seq_len(3)), observer = factor(seq_len(2))
  )
  block_effect <- c(-0.8, 0.1, 0.7)
  observer_effect <- c(-0.4, 0.4)
  dat$y <- 1.5 * (dat$treatment == "impact") +
    block_effect[dat$block] + observer_effect[dat$observer] +
    rep(c(-0.1, 0.1), length.out = nrow(dat))
  dat
}

.blocked_permutation_example <- function(treatments = 2L, blocks = 3L,
                                         replicates = 2L, seed = 91L) {
  dat <- expand.grid(
    replicate = factor(seq_len(replicates)),
    treatment = factor(seq_len(treatments)),
    block = factor(seq_len(blocks))
  )
  set.seed(seed)
  cell <- interaction(dat$treatment, dat$block, drop = TRUE)
  dat$y <- rnorm(blocks, sd = 0.7)[dat$block] +
    rnorm(treatments * blocks, sd = 1.2)[cell] +
    rnorm(nrow(dat), sd = 0.4)
  dat
}

.direct_blocked_f <- function(response, index, dat) {
  reference <- summary(aov(
    response[index] ~ treatment * block, dat
  ))[[1]]
  rows <- trimws(rownames(reference))
  reference[["Mean Sq"]][match("treatment", rows)] /
    reference[["Mean Sq"]][match("treatment:block", rows)]
}

test_that("crossed mixed engine derives exact EMS denominator and components", {
  dat <- .crossed_example()
  design <- compile_permanova_design(
    ~ treatment + (1 | block) + (1 | observer), dat, distance = dist(dat$y)
  )
  expect_error(permanova_crossed(design), "Unmodelled fixed-by-random")
  fit <- permanova_crossed(design, truncate = TRUE, assume_additive = TRUE)

  expect_s3_class(fit, "permanova_crossed")
  expect_s3_class(fit$tests$treatment, "permanova_test_spec")
  expect_true(is.na(fit$tests$treatment$p_value))
  expect_match(fit$tests$treatment$approximation_status, "exact EMS denominator")
  expect_true(fit$assumptions$additive_opt_in)
  expect_equal(sort(fit$assumptions$zero_variance_components),
               c("block:treatment", "observer:treatment"))
  expect_match(fit$tests$treatment$symbolic_ems$derivation, "zero variance")
  expect_gt(length(fit$tests$treatment$denominator$pooled_strata), 1)
  expect_named(fit$variance_components$unrestricted, c("residual", "block", "observer"))
  expect_true(all(fit$variance_components$truncated >= 0))
  expect_equal(
    as.vector(fit$variance_components$coefficient_matrix %*%
      fit$variance_components$unrestricted),
    unname(fit$variance_components$observed_mean_squares), tolerance = 1e-10
  )

  reference <- summary(aov(y ~ treatment + block + observer, dat))[[1]]
  expect_equal(
    fit$tests$treatment$statistic,
    unname(reference[["F value"]][[1L]]), tolerance = 1e-10
  )
})

test_that("Satterthwaite p-values require explicit opt-in", {
  dat <- .crossed_example()
  design <- compile_permanova_design(
    ~ treatment + (1 | block) + (1 | observer), dat, distance = dist(dat$y)
  )
  fit <- permanova_crossed(
    design, p_method = "satterthwaite", assume_additive = TRUE
  )
  expect_true(is.finite(fit$tests$treatment$p_value))
  expect_true(any(grepl("approximate", fit$tests$treatment$warnings)))
})

test_that("negative unrestricted variance components are retained", {
  dat <- .crossed_example()
  dat$y <- ifelse(dat$replicate == 1, -1, 1)
  design <- compile_permanova_design(
    ~ treatment + (1 | block) + (1 | observer), dat, distance = dist(dat$y)
  )
  components <- estimate_variance_components(
    design, truncate = TRUE, assume_additive = TRUE
  )
  expect_equal(components$truncated, pmax(components$unrestricted, 0))
  expect_true(any(components$unrestricted < 0))
  expect_true(any(grepl("retained", components$warnings)))
})

test_that("unbalanced crossed designs fail explicitly", {
  dat <- .crossed_example()[-1, ]
  design <- compile_permanova_design(
    ~ treatment + (1 | block) + (1 | observer), dat, distance = dist(dat$y)
  )
  expect_error(permanova_crossed(design), "complete, equally replicated")
})

test_that("declared fixed-by-random interaction supplies the exact denominator", {
  dat <- expand.grid(
    replicate = factor(1:2), treatment = factor(c("control", "impact")),
    block = factor(1:6)
  )
  set.seed(12)
  block_effect <- rnorm(6)
  interaction_effect <- rnorm(12, sd = 2)
  cell <- interaction(dat$treatment, dat$block, drop = TRUE)
  dat$y <- block_effect[dat$block] + interaction_effect[cell] + rnorm(nrow(dat))
  design <- compile_permanova_design(
    ~ treatment + (1 | block) + (1 | treatment:block), dat, dist(dat$y)
  )
  fit <- permanova_crossed(design, p_method = "satterthwaite")
  reference <- summary(aov(y ~ treatment * block, dat))[[1]]
  rows <- trimws(rownames(reference))
  classical_f <- reference[["Mean Sq"]][match("treatment", rows)] /
    reference[["Mean Sq"]][match("treatment:block", rows)]

  expect_equal(fit$tests$treatment$statistic, unname(classical_f), tolerance = 1e-10)
  expect_equal(fit$tests$treatment$denominator$pooled_strata, "treatment:block")
  expect_equal(fit$tests$treatment$denominator$expression,
               "1*MS[treatment:block]")
  expect_equal(fit$tests$treatment$symbolic_ems$zero_variance_assumptions,
               character())
  expect_false(fit$assumptions$additive_opt_in)
})

test_that("non-Euclidean Satterthwaite p-values are suppressed", {
  dat <- expand.grid(
    replicate = factor(1:2), treatment = factor(c("a", "b")), block = factor(1:4)
  )
  set.seed(8)
  points <- matrix(rnorm(nrow(dat) * 2), ncol = 2)
  D <- as.matrix(dist(points))
  D[1, 2] <- D[2, 1] <- 4 * max(D)
  design <- compile_permanova_design(
    ~ treatment + (1 | block) + (1 | treatment:block), dat, D
  )
  expect_false(design$spectral$euclidean)
  expect_warning(
    fit <- permanova_crossed(design, p_method = "satterthwaite"),
    "suppressed"
  )
  expect_true(is.na(fit$tests$treatment$p_value))
  expect_match(fit$tests$treatment$approximation_status, "suppressed")
  expect_match(fit$tests$treatment$warnings, "non-Euclidean")
})

test_that("explicit crossed random interactions reach EMS inference", {
  dat <- expand.grid(
    replicate = factor(1:2), treatment = factor(c("a", "b")),
    block = factor(1:3), observer = factor(1:2)
  )
  set.seed(21)
  dat$y <- rnorm(nrow(dat))
  design <- compile_permanova_design(
    ~ treatment +
      (1 | block) + (1 | observer) + (1 | block:observer) +
      (1 | treatment:block) + (1 | treatment:observer) +
      (1 | treatment:block:observer),
    dat, dist(dat$y)
  )
  fit <- permanova_crossed(design)

  expect_s3_class(fit, "permanova_crossed")
  expect_equal(nrow(design$containment_graph), 0)
  expect_equal(design$random_term_table$origin[[3L]], "interaction")
  expect_true("block:observer" %in% names(fit$variance_components$unrestricted))
  expect_true("treatment:block:observer" %in%
                fit$tests$treatment$denominator$pooled_strata)
})

test_that("two-treatment crossed permutations enumerate cells within blocks", {
  dat <- .blocked_permutation_example(treatments = 2, blocks = 3)
  design <- compile_permanova_design(
    ~ treatment + (1 | block) + (1 | treatment:block),
    dat, dist(dat$y)
  )
  fit <- suppressWarnings(permanova_crossed(
    design, p_method = "permutation", permutations = "all",
    keep_permutations = TRUE
  ))

  expect_equal(fit$permutation$attainable, 2 ^ 3)
  expect_equal(fit$permutation$used, 2 ^ 3)
  expect_true(fit$permutation$exact)
  expect_equal(fit$tests$treatment$unique_permutations_used, 2 ^ 3)
  expect_equal(fit$tests$treatment$exchangeable_units, "treatment:block")
  expect_equal(fit$tests$treatment$permutation_strata, "block")
  expect_equal(fit$tests$treatment$p_value * 2 ^ 3 %% 1,
               0, tolerance = 1e-12)
  expect_match(fit$tests$treatment$approximation_status, "exhaustive")

  structure <- distmer:::.validate_crossed_design(design)
  permutation_structure <-
    distmer:::.validate_crossed_permutation_design(design, structure)
  plan <- distmer:::.crossed_block_permutation_plan(
    permutation_structure, "all", 100L
  )
  direct <- vapply(seq_len(nrow(plan$orders)), function(row) {
    index <- distmer:::.block_permutation_index(
      permutation_structure$cell_blocks, plan$orders[row, ]
    )
    .direct_blocked_f(dat$y, index, dat)
  }, numeric(1))
  expect_equal(fit$tests$treatment$permutation_statistics, direct,
               tolerance = 1e-10)
  for (cells in permutation_structure$cells_by_block) {
    expect_true(all(apply(plan$orders[, cells, drop = FALSE], 1, function(x) {
      identical(sort(as.integer(x)), sort(as.integer(cells)))
    })))
  }
  for (row in seq_len(nrow(plan$orders))) {
    index <- distmer:::.block_permutation_index(
      permutation_structure$cell_blocks, plan$orders[row, ]
    )
    for (target in seq_along(permutation_structure$cell_blocks)) {
      expect_identical(
        index[permutation_structure$cell_blocks[[target]]],
        permutation_structure$cell_blocks[[plan$orders[row, target]]]
      )
    }
  }
})

test_that("three-treatment crossed permutations enumerate the full product", {
  dat <- .blocked_permutation_example(treatments = 3, blocks = 2)
  design <- compile_permanova_design(
    ~ treatment + (1 | block) + (1 | block:treatment),
    dat, dist(dat$y)
  )
  fit <- suppressWarnings(permanova_crossed(
    design, p_method = "permutation", permutations = "all",
    keep_permutations = TRUE
  ))
  expect_equal(fit$permutation$attainable, factorial(3) ^ 2)
  expect_equal(fit$tests$treatment$unique_permutations_used, factorial(3) ^ 2)

  structure <- distmer:::.validate_crossed_design(design)
  permutation_structure <-
    distmer:::.validate_crossed_permutation_design(design, structure)
  plan <- distmer:::.crossed_block_permutation_plan(
    permutation_structure, "all", 100L
  )
  direct <- vapply(seq_len(nrow(plan$orders)), function(row) {
    index <- distmer:::.block_permutation_index(
      permutation_structure$cell_blocks, plan$orders[row, ]
    )
    .direct_blocked_f(dat$y, index, dat)
  }, numeric(1))
  expect_equal(fit$tests$treatment$permutation_statistics, direct,
               tolerance = 1e-10)
})

test_that("crossed permutation p-values retain non-Euclidean spectra", {
  dat <- .blocked_permutation_example()
  set.seed(93)
  D <- as.matrix(dist(matrix(rnorm(nrow(dat) * 2), ncol = 2)))
  D[1, 2] <- D[2, 1] <- 4 * max(D)
  design <- compile_permanova_design(
    ~ treatment + (1 | block) + (1 | treatment:block), dat, D
  )
  expect_false(design$spectral$euclidean)
  fit <- suppressWarnings(permanova_crossed(
    design, p_method = "permutation", permutations = "all"
  ))
  expect_true(is.finite(fit$tests$treatment$p_value))
  expect_gte(fit$tests$treatment$p_value, 0)
  expect_lte(fit$tests$treatment$p_value, 1)
  expect_true(any(grepl(
    "retained non-Euclidean", fit$tests$treatment$warnings
  )))
  expect_equal(fit$design$spectral$negative_count,
               design$spectral$negative_count)
})

test_that("crossed Monte Carlo permutations are seeded and validated", {
  dat <- .blocked_permutation_example(blocks = 8)
  design <- compile_permanova_design(
    ~ treatment + (1 | block) + (1 | treatment:block),
    dat, dist(dat$y)
  )
  first <- permanova_crossed(
    design, p_method = "permutation", permutations = 19,
    seed = 44, keep_permutations = TRUE
  )
  second <- permanova_crossed(
    design, p_method = "permutation", permutations = 19,
    seed = 44, keep_permutations = TRUE
  )
  expect_false(first$permutation$exact)
  expect_equal(first$permutation$attainable, 2 ^ 8)
  expect_equal(first$tests$treatment$permutation_statistics,
               second$tests$treatment$permutation_statistics)
  expect_equal(first$tests$treatment$p_value,
               second$tests$treatment$p_value)
  expect_equal(first$tests$treatment$p_value * 20 %% 1,
               0, tolerance = 1e-12)
  expect_error(permanova_crossed(
    design, p_method = "permutation", permutations = -1
  ), "non-negative integer")
  expect_error(permanova_crossed(
    design, p_method = "permutation", permutations = "all", max_exact = 100
  ), "exceeds `max_exact`")
})

test_that("crossed permutation envelope rejects unsupported designs", {
  dat <- .blocked_permutation_example()
  missing_interaction <- compile_permanova_design(
    ~ treatment + (1 | block), dat, dist(dat$y)
  )
  expect_error(permanova_crossed(
    missing_interaction, p_method = "permutation"
  ), "Unmodelled fixed-by-random")
  expect_error(permanova_crossed(
    missing_interaction, p_method = "permutation", assume_additive = TRUE
  ), "does not permit")

  nested <- compile_permanova_design(
    ~ treatment + (1 | block/treatment), dat, dist(dat$y)
  )
  expect_error(permanova_crossed(
    nested, p_method = "permutation"
  ), "/`-generated nesting")

  extra <- expand.grid(
    replicate = factor(1:2), treatment = factor(1:2),
    block = factor(1:3), observer = factor(1:2)
  )
  set.seed(94)
  extra$y <- rnorm(nrow(extra))
  extra_design <- compile_permanova_design(
    ~ treatment + (1 | block) + (1 | observer) +
      (1 | treatment:block) + (1 | treatment:observer),
    extra, dist(extra$y)
  )
  expect_error(permanova_crossed(
    extra_design, p_method = "permutation"
  ), "exactly one random block")

  one_replicate <- subset(dat, replicate == levels(replicate)[[1L]])
  one_replicate$replicate <- droplevels(one_replicate$replicate)
  one_design <- compile_permanova_design(
    ~ treatment + (1 | block) + (1 | treatment:block),
    one_replicate, dist(one_replicate$y)
  )
  expect_error(permanova_crossed(
    one_design, p_method = "permutation"
  ), "At least two observations")

  unbalanced <- dat[-1L, ]
  unbalanced_design <- compile_permanova_design(
    ~ treatment + (1 | block) + (1 | treatment:block),
    unbalanced, dist(unbalanced$y)
  )
  expect_error(permanova_crossed(
    unbalanced_design, p_method = "permutation"
  ), "complete, equally replicated")

  continuous <- dat
  continuous$treatment_numeric <- as.numeric(continuous$treatment)
  continuous_design <- compile_permanova_design(
    ~ treatment_numeric + (1 | block) +
      (1 | treatment_numeric:block),
    continuous, dist(continuous$y)
  )
  expect_error(permanova_crossed(
    continuous_design, p_method = "permutation"
  ), "categorical fixed and random factors")

  multifactor <- expand.grid(
    replicate = factor(1:2), treatment = factor(1:2),
    habitat = factor(1:2), block = factor(1:3)
  )
  set.seed(95)
  multifactor$y <- rnorm(nrow(multifactor))
  multifactor_design <- compile_permanova_design(
    ~ treatment + habitat + (1 | block) +
      (1 | treatment:block) + (1 | habitat:block) +
      (1 | treatment:habitat:block),
    multifactor, dist(multifactor$y)
  )
  expect_error(permanova_crossed(
    multifactor_design, p_method = "permutation"
  ), "exactly one categorical fixed treatment")
})
