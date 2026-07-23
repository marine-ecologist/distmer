.phase19_reference_f <- function(y, treatment, contrast, block = NULL) {
  data <- data.frame(treatment = treatment)
  formula <- ~ treatment
  if (!is.null(block)) {
    data$block <- block
    formula <- ~ block + treatment
  }
  X <- model.matrix(formula, data)
  beta <- qr.solve(X, y)
  treatment_levels <- levels(treatment)
  if (is.null(block)) {
    reference <- data.frame(treatment = factor(
      treatment_levels, levels = treatment_levels
    ))
  } else {
    reference <- expand.grid(
      treatment = treatment_levels, block = levels(block),
      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
    )
    reference$treatment <- factor(reference$treatment,
                                  levels = treatment_levels)
    reference$block <- factor(reference$block, levels = levels(block))
  }
  reference_X <- model.matrix(formula, reference)
  means <- do.call(rbind, lapply(treatment_levels, function(level) {
    colMeans(reference_X[reference$treatment == level, , drop = FALSE])
  }))
  L <- as.numeric(contrast %*% means)
  estimate <- sum(L * beta)
  numerator <- estimate ^ 2 / as.numeric(L %*% solve(crossprod(X), L))
  residual <- y - X %*% beta
  denominator <- sum(residual ^ 2) / (length(y) - qr(X)$rank)
  numerator / denominator
}

.phase19_unequal_random <- function() {
  sizes <- rep(2:4, each = 3)
  treatment_by_site <- rep(c("A", "B", "C"), 3)
  data <- data.frame(
    site = factor(rep(seq_along(sizes), sizes)),
    treatment = factor(rep(treatment_by_site, sizes),
                       levels = c("A", "B", "C"))
  )
  site_effect <- c(-1, 0.2, 1, -0.4, 0.7, 1.4, 0.1, 1.3, 2.8)
  data$y <- site_effect[data$site] +
    ave(seq_len(nrow(data)), data$site, FUN = function(index) {
      seq(-0.15, 0.15, length.out = length(index))
    })
  data
}

.phase19_nested <- function() {
  topology <- list(
    c(2L, 3L), c(2L, 3L), c(2L, 3L),
    c(2L, 2L, 2L), c(2L, 2L, 2L), c(2L, 2L, 2L)
  )
  treatment <- rep(c("A", "B", "C"), 2)
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
  data$treatment <- factor(data$treatment, levels = c("A", "B", "C"))
  region_effect <- c(-1, 0, 1, -0.2, 1, 2.4)
  data$y <- region_effect[data$region] + as.integer(data$site) / 30 +
    data$replicate / 100
  data
}

.phase19_block <- function() {
  cells <- data.frame(
    block = factor(c(1, 1, 1, 2, 2, 2, 3, 3, 4, 4)),
    treatment = factor(
      c("A", "B", "C", "A", "B", "C", "A", "B", "B", "C"),
      levels = c("A", "B", "C")
    ),
    size = c(2L, 2L, 2L, 3L, 3L, 3L, 2L, 2L, 2L, 2L)
  )
  data <- cells[rep(seq_len(nrow(cells)), cells$size),
                c("block", "treatment")]
  rownames(data) <- NULL
  data$replicate <- ave(
    seq_len(nrow(data)), interaction(data$block, data$treatment),
    FUN = seq_along
  )
  cell <- interaction(data$treatment, data$block, drop = TRUE)
  cell_effect <- c(-.7, .2, 1.1, .4, -.3, .8, .1, -.5, .9, 1.5)
  data$y <- c(-.5, .2, .7, 1.1)[data$block] +
    cell_effect[cell] + data$replicate / 100
  data
}

test_that("unequal random-intercept contrasts match unit-centroid GLH", {
  data <- .phase19_unequal_random()
  design <- compile_permanova_design(
    ~ treatment + (1 | site), data, dist(data$y)
  )
  contrast <- c(A = 1, B = -1, C = 0)
  fit <- permanova_mixed_contrasts(design, contrast, permutations = 0)
  means <- aggregate(y ~ site + treatment, data, mean)
  reference <- .phase19_reference_f(
    means$y, means$treatment, contrast
  )

  expect_equal(fit$table$pseudo_f, reference, tolerance = 1e-10)
  expect_match(fit$design_type, "unequally replicated")
  expect_match(fit$tests[[1]]$permutation_strata, "observation count")
})

test_that("unbalanced nested contrasts match top-centroid GLH", {
  data <- .phase19_nested()
  design <- compile_permanova_design(
    ~ treatment + (1 | region/site), data, dist(data$y)
  )
  contrast <- c(A = 1, B = -1, C = 0)
  fit <- permanova_mixed_contrasts(design, contrast, permutations = 0)
  means <- aggregate(y ~ region + treatment, data, mean)
  reference <- .phase19_reference_f(
    means$y, means$treatment, contrast
  )

  expect_equal(fit$table$pseudo_f, reference, tolerance = 1e-10)
  expect_match(fit$design_type, "unbalanced top-level")
  expect_match(fit$tests[[1]]$permutation_strata, "child observation-count")
})

test_that("incomplete-block contrasts use common-grid adjusted means", {
  data <- .phase19_block()
  design <- compile_permanova_design(
    ~ treatment + (1 | block) + (1 | treatment:block),
    data, dist(data$y)
  )
  contrast <- c(A = 1, B = -1, C = 0)
  fit <- permanova_mixed_contrasts(design, contrast, permutations = 0)
  structure <- distmer:::.validate_unbalanced_block(design)
  means <- aggregate(y ~ treatment + block, data, mean)
  means <- means[match(
    paste(structure$cell_mapping$treatment, structure$cell_mapping$block),
    paste(means$treatment, means$block)
  ), ]
  reference <- .phase19_reference_f(
    means$y, means$treatment, contrast, means$block
  )

  expect_equal(fit$table$pseudo_f, reference, tolerance = 1e-10)
  expect_match(fit$design_type, "incomplete or unequally")
  expect_match(fit$estimand, "additive block adjustment")
})

test_that("unbalanced mixed contrast permutations equal direct reconstruction", {
  data <- .phase19_block()
  design <- compile_permanova_design(
    ~ treatment + (1 | block) + (1 | treatment:block),
    data, dist(data$y)
  )
  contrast <- c(A = 1, B = -1, C = 0)
  fit <- permanova_mixed_contrasts(
    design, contrast, permutations = "all", keep_permutations = TRUE
  )
  structure <- distmer:::.mixed_contrast_structure(design)
  units <- distmer:::.prepare_mixed_contrast_units(structure)
  plan <- distmer:::.mixed_contrast_plan(structure, "all", 100000)
  coefficient_row <- as.numeric(contrast %*% units$mean_rows)
  projector <- distmer:::.linear_contrast_projector(
    units$X, coefficient_row, units$H_full, 1e-10
  )
  y <- vapply(structure$aggregated$blocks, function(index) mean(data$y[index]),
              numeric(1))
  residual <- diag(length(y)) - projector$H_reduced
  direct <- apply(plan$orders, 1, function(order) {
    permuted <- (projector$H_reduced + residual[order, , drop = FALSE]) %*% y
    sum(permuted * (projector$H %*% permuted)) /
      (sum(permuted * (units$H_denominator %*% permuted)) /
         units$denominator_df)
  })

  expect_equal(fit$tests[[1]]$permutation_statistics, direct,
               tolerance = 1e-10)
  expect_true(fit$permutation$exact)
})

test_that("unbalanced mixed pairwise families retain multiplicity adjustment", {
  data <- .phase19_unequal_random()
  design <- compile_permanova_design(
    ~ treatment + (1 | site), data, dist(data$y)
  )
  fit <- permanova_mixed_pairwise(
    design, permutations = 19, seed = 5, adjust = "holm"
  )

  expect_equal(nrow(fit$table), 3)
  expect_equal(fit$table$p_adjusted, p.adjust(fit$table$p_perm, "holm"))
})
