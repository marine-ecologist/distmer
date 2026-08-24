.unbalanced_nested_example <- function() {
  topology <- list(c(2L, 3L), c(2L, 3L), c(2L, 2L, 2L),
                   c(2L, 2L, 2L), c(3L, 4L), c(3L, 4L))
  treatment <- rep(c("control", "impact"), 3L)
  rows <- list()
  child_id <- 0L
  for (parent in seq_along(topology)) {
    for (size in topology[[parent]]) {
      child_id <- child_id + 1L
      rows[[length(rows) + 1L]] <- data.frame(
        region = factor(parent, levels = seq_along(topology)),
        site = factor(child_id), treatment = treatment[[parent]],
        replicate = seq_len(size)
      )
    }
  }
  dat <- do.call(rbind, rows)
  dat$region <- factor(dat$region)
  dat$site <- factor(dat$site)
  dat$treatment <- factor(dat$treatment, levels = c("control", "impact"))
  region_effect <- c(-0.9, 0.2, -0.4, 0.8, -0.2, 1.1)
  site_effect <- seq(-0.35, 0.35, length.out = nlevels(dat$site))
  dat$y <- 0.5 * (dat$treatment == "impact") +
    region_effect[dat$region] + site_effect[dat$site] +
    ave(dat$replicate, dat$site, FUN = function(x) {
      seq(-0.1, 0.1, length.out = length(x))
    })
  dat
}

.unbalanced_block_example <- function() {
  cells <- data.frame(
    block = factor(c(1, 1, 1, 2, 2, 2, 3, 3, 4, 4)),
    treatment = factor(c("A", "B", "C", "A", "B", "C", "A", "B", "B", "C"),
                       levels = c("A", "B", "C")),
    size = c(2L, 2L, 2L, 3L, 3L, 3L, 2L, 2L, 2L, 2L)
  )
  dat <- cells[rep(seq_len(nrow(cells)), cells$size), c("block", "treatment")]
  rownames(dat) <- NULL
  dat$replicate <- ave(seq_len(nrow(dat)), interaction(dat$block, dat$treatment),
                       FUN = seq_along)
  block_effect <- c(-0.6, 0.1, 0.4, -0.2)
  cell_effect <- c(-0.5, 0.2, 0.9, 0.7, -0.4, 0.3, 0.1, -0.8, 0.6, -0.2)
  cell_key <- paste(cells$block, cells$treatment)
  observation_cell <- match(paste(dat$block, dat$treatment), cell_key)
  cell <- interaction(dat$treatment, dat$block, drop = TRUE)
  dat$y <- c(A = 0, B = 0.4, C = 0.8)[dat$treatment] +
    block_effect[dat$block] + cell_effect[observation_cell] +
    ave(dat$replicate, cell, FUN = function(x) {
      seq(-0.12, 0.12, length.out = length(x))
    })
  dat
}

test_that("unbalanced nested F and exhaustive permutations match top means", {
  dat <- .unbalanced_nested_example()
  design <- compile_permanova_design(
    ~ treatment + (1 | region/site), dat, dist(dat$y)
  )
  fit <- suppressWarnings(permanova_unbalanced_nested(
    design, permutations = "all", keep_permutations = TRUE
  ))
  # parent centroids average their child centroids with equal weight
  site_means <- aggregate(y ~ site + region + treatment, dat, mean)
  means <- aggregate(y ~ region + treatment, site_means, mean)
  means <- means[order(as.integer(means$region)), ]
  reference <- summary(aov(y ~ treatment, means))[[1]][1L, "F value"]
  expect_equal(fit$test$statistic, unname(reference), tolerance = 1e-10)
  expect_equal(fit$permutation$attainable, 8)

  structure <- .validate_unbalanced_nested(design)
  plan <- .stratified_permutation_plan(structure$topology_strata, "all", 100L)
  reduced <- lm(y ~ 1, means)
  direct <- apply(plan$orders, 1L, function(index) {
    response <- fitted(reduced) + residuals(reduced)[index]
    summary(aov(response ~ treatment, means))[[1]][1L, "F value"]
  })
  expect_equal(fit$test$permutation_statistics, direct, tolerance = 1e-10)
})

test_that("unbalanced block F and cell permutations match scalar references", {
  dat <- .unbalanced_block_example()
  design <- compile_permanova_design(
    ~ treatment + (1 | block) + (1 | treatment:block), dat, dist(dat$y)
  )
  fit <- suppressWarnings(permanova_unbalanced_block(
    design, permutations = "all", keep_permutations = TRUE
  ))
  means <- aggregate(y ~ treatment + block, dat, mean)
  structure <- .validate_unbalanced_block(design)
  means <- means[match(
    paste(structure$cell_mapping$treatment, structure$cell_mapping$block),
    paste(means$treatment, means$block)
  ), ]
  reference <- summary(aov(y ~ block + treatment, means))[[1]]
  reference_f <- reference["treatment", "Mean Sq"] / reference["Residuals", "Mean Sq"]
  expect_equal(fit$test$statistic, unname(reference_f), tolerance = 1e-10)
  expect_equal(fit$permutation$attainable, 144)

  plan <- .stratified_permutation_plan(
    structure$exchangeability_strata, "all", 1000L
  )
  direct <- apply(plan$orders, 1L, function(index) {
    response <- means$y[index]
    table <- summary(aov(response ~ block + treatment, means))[[1]]
    table["treatment", "Mean Sq"] / table["Residuals", "Mean Sq"]
  })
  expect_equal(fit$test$permutation_statistics, direct, tolerance = 1e-10)
})

test_that("Phase 12 rejects topology and cell-size confounding", {
  nested <- .unbalanced_nested_example()
  nested$treatment <- factor(ifelse(
    nested$region %in% c("1", "2"), "control", "impact"
  ))
  nested_design <- compile_permanova_design(
    ~ treatment + (1 | region/site), nested, dist(nested$y)
  )
  expect_error(
    permanova_unbalanced_nested(nested_design, permutations = 0),
    "do not span every treatment contrast"
  )

  blocked <- .unbalanced_block_example()
  omitted <- compile_permanova_design(
    ~ treatment + (1 | block), blocked, dist(blocked$y)
  )
  expect_error(permanova_unbalanced_block(omitted, permutations = 0),
               "explicit treatment:block")

  cell <- interaction(blocked$treatment, blocked$block, drop = TRUE)
  keep <- ave(seq_len(nrow(blocked)), cell, FUN = seq_along) <=
    as.integer(blocked$treatment)
  confounded <- blocked[keep, ]
  confounded_design <- compile_permanova_design(
    ~ treatment + (1 | block) + (1 | treatment:block),
    confounded, dist(confounded$y)
  )
  expect_error(
    permanova_unbalanced_block(confounded_design, permutations = 0),
    "do not span every treatment contrast"
  )
})

test_that("Phase 12 rejects fully topology-confounded nested layouts", {
  topology <- list(
    c(2L, 2L), c(2L, 2L), c(2L, 2L),
    c(3L, 3L), c(3L, 3L), c(3L, 3L)
  )
  treatment <- c("control", "control", "control",
                 "impact", "impact", "impact")
  rows <- list()
  site <- 0L
  for (region in seq_along(topology)) {
    for (size in topology[[region]]) {
      site <- site + 1L
      rows[[length(rows) + 1L]] <- data.frame(
        region = region, site = site, treatment = treatment[[region]],
        replicate = seq_len(size)
      )
    }
  }
  data <- do.call(rbind, rows)
  data$region <- factor(data$region)
  data$site <- factor(data$site)
  data$treatment <- factor(data$treatment)
  data$y <- as.integer(data$region) + as.integer(data$site) / 20 +
    data$replicate / 100
  design <- compile_permanova_design(
    ~ treatment + (1 | region/site), data, dist(data$y)
  )

  expect_error(
    permanova_unbalanced_nested(design, permutations = 0),
    "do not span every treatment contrast"
  )
})

test_that("Phase 12 rejects treatment-confounded within-block size strata", {
  cells <- expand.grid(
    treatment = factor(c("A", "B", "C"), levels = c("A", "B", "C")),
    block = factor(seq_len(4L))
  )
  cells$size <- c(A = 2L, B = 3L, C = 4L)[cells$treatment]
  data <- cells[rep(seq_len(nrow(cells)), cells$size),
                c("treatment", "block")]
  rownames(data) <- NULL
  data$replicate <- ave(seq_len(nrow(data)),
                        interaction(data$treatment, data$block),
                        FUN = seq_along)
  data$y <- as.integer(data$block) +
    c(A = 0, B = 0.5, C = 1)[data$treatment] + data$replicate / 100
  design <- compile_permanova_design(
    ~ treatment + (1 | block) + (1 | treatment:block),
    data, dist(data$y)
  )

  expect_error(
    permanova_unbalanced_block(design, permutations = 0),
    "do not span every treatment contrast"
  )
})

test_that("Phase 12 retains non-Euclidean empirical p-values", {
  dat <- .unbalanced_block_example()
  D <- as.matrix(dist(cbind(dat$y, seq_len(nrow(dat)) / 20)))
  D[1L, 2L] <- D[2L, 1L] <- 5 * max(D)
  design <- compile_permanova_design(
    ~ treatment + (1 | block) + (1 | treatment:block), dat, D
  )
  fit <- suppressWarnings(permanova_unbalanced_block(
    design, permutations = 19
  ))
  expect_false(design$spectral$euclidean)
  expect_true(is.finite(fit$test$p_value))
  expect_true(any(grepl("non-Euclidean", fit$test$warnings)))
})

test_that("a permutation stratum holding one unit is kept fixed, not resampled", {
  # `sample(x)` permutes seq_len(x) when x is a single number, so a singleton
  # stratum used to crash the restricted plan (or silently misbehave when the
  # lone index happened to be 1). Every engine using restricted permutation
  # depends on this.
  strata <- factor(c("a", "a", "b", "c", "c", "c"))
  positions <- unname(split(seq_along(strata), strata))
  set.seed(1)
  plan <- .stratified_permutation_plan(strata, 50L, 100000L)

  expect_gt(nrow(plan$orders), 0L)
  for (row in seq_len(nrow(plan$orders))) {
    for (block in positions) {
      # a unit may move only inside its own stratum
      expect_setequal(plan$orders[row, block], block)
    }
    # the singleton stratum can only map to itself
    expect_identical(plan$orders[row, positions[[2L]]], positions[[2L]])
  }

  # and the singleton index must not be 1, or the bug would be masked
  expect_gt(positions[[2L]][[1L]], 1L)
})

test_that("an engine runs when a unit-size stratum is a singleton", {
  set.seed(4)
  rows <- list()
  sizes <- c(2L, 2L, 3L, 3L, 4L, 5L)   # sizes 4 and 5 are singleton strata
  for (cluster in seq_along(sizes)) {
    rows[[length(rows) + 1L]] <- data.frame(
      unit = paste0("u", cluster), rep = seq_len(sizes[[cluster]]),
      stringsAsFactors = FALSE
    )
  }
  d <- do.call(rbind, rows)
  d$unit <- factor(d$unit)
  d$tr <- factor(ifelse(as.integer(d$unit) %% 2L == 0L, "control", "impact"))
  d$y <- stats::rnorm(nrow(d))

  expect_no_error(
    permanova_unequal_random_intercept(
      compile_permanova_design(~ tr + (1 | unit), d, dist(d$y)),
      permutations = 49L, seed = 1L
    )
  )
})
