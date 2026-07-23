.lower_child_example <- function(zero_residual = FALSE) {
  rows <- list()
  child_id <- 0L
  for (parent in 1:2) {
    for (size in 2:3) {
      for (treatment in c("A", "B")) {
        child_id <- child_id + 1L
        rows[[length(rows) + 1L]] <- data.frame(
          region = parent, site = child_id, treatment = treatment,
          observation = seq_len(size)
        )
      }
    }
  }
  data <- do.call(rbind, rows)
  data$region <- factor(data$region)
  data$site <- factor(data$site)
  data$treatment <- factor(data$treatment, levels = c("A", "B"))
  child_effect <- if (zero_residual) rep(0, nlevels(data$site)) else
    seq(-0.35, 0.35, length.out = nlevels(data$site))
  data$y <- c(`1` = -0.4, `2` = 0.5)[data$region] +
    c(A = -0.3, B = 0.6)[data$treatment] +
    child_effect[data$site] +
    if (zero_residual) 0 else sin(seq_len(nrow(data))) / 20
  data
}

.lower_observation_example <- function() {
  parent_by_child <- c(1, 1, 2, 2, 2)
  data <- do.call(rbind, lapply(seq_along(parent_by_child), function(child) {
    data.frame(
      region = parent_by_child[[child]], site = child,
      treatment = c("A", "B"), observation = 1:2
    )
  }))
  data$region <- factor(data$region)
  data$site <- factor(data$site)
  data$treatment <- factor(data$treatment, levels = c("A", "B"))
  data$y <- seq(-0.3, 0.4, length.out = nlevels(data$site))[data$site] +
    c(A = -0.5, B = 0.7)[data$treatment] +
    cos(seq_len(nrow(data))) / 20
  data
}

.lower_design <- function(data) {
  compile_permanova_design(
    ~ treatment + (1 | region / site), data, dist(data$y)
  )
}

test_that("child-level unbalanced nested F equals child-centroid GLH", {
  data <- .lower_child_example()
  design <- .lower_design(data)
  fit <- permanova_unbalanced_nested_lower(design, permutations = 0)
  centroids <- aggregate(y ~ site + region + treatment, data, mean)
  reference <- anova(lm(y ~ region + treatment, centroids))["treatment", "F value"]

  expect_equal(fit$assignment_level, "child")
  expect_equal(fit$table$pseudo_f[[1]], reference, tolerance = 1e-10)
  expect_match(fit$estimand, "equally weighted.*site.*parent|equally weighted.*site.*region")
  expect_match(fit$test$permutation_strata, "exact.*observation count")
})

test_that("observation-level unbalanced nested F equals child-adjusted GLH", {
  data <- .lower_observation_example()
  design <- .lower_design(data)
  fit <- permanova_unbalanced_nested_lower(design, permutations = 0)
  reference <- anova(lm(y ~ site + treatment, data))["treatment", "F value"]

  expect_equal(fit$assignment_level, "observation")
  expect_equal(fit$table$pseudo_f[[1]], reference, tolerance = 1e-10)
  expect_match(fit$estimand, "observation-level")
  expect_equal(fit$test$permutation_strata, "region:site")
})

test_that("lower-level restricted permutation vectors equal reconstruction", {
  for (data in list(.lower_child_example(), .lower_observation_example())) {
    design <- .lower_design(data)
    fit <- permanova_unbalanced_nested_lower(
      design, permutations = "all", keep_permutations = TRUE,
      max_exact = 1000
    )
    validated <- distmer:::.validate_unbalanced_nested_lower(design)
    plan <- distmer:::.mixed_contrast_plan(
      validated$structure, "all", 1000
    )
    y <- vapply(
      validated$structure$aggregated$blocks,
      function(index) mean(data$y[index]), numeric(1)
    )
    residual <- diag(length(y)) - validated$H_reduced
    direct <- apply(plan$orders, 1L, function(order) {
      permuted <- (validated$H_reduced +
        residual[order, , drop = FALSE]) %*% y
      sum(permuted * (validated$H_numerator %*% permuted)) /
        validated$numerator_df /
        (sum(permuted * (validated$units$H_denominator %*% permuted)) /
           validated$units$denominator_df)
    })
    expect_true(fit$permutation$exact)
    expect_equal(fit$test$permutation_statistics, direct,
                 tolerance = 1e-10)
  }
})

test_that("lower-level mixed contrasts reuse the same unit GLH", {
  for (data in list(.lower_child_example(), .lower_observation_example())) {
    design <- .lower_design(data)
    omnibus <- permanova_unbalanced_nested_lower(design, permutations = 0)
    contrast <- permanova_mixed_contrasts(
      design, c(A = 1, B = -1), permutations = 0
    )
    pairwise <- permanova_mixed_pairwise(
      design, permutations = 15, seed = 3
    )
    maximum <- permanova_mixed_maxT(
      design, c(A = 1, B = -1), permutations = 15, seed = 3
    )

    expect_equal(contrast$table$pseudo_f, omnibus$table$pseudo_f[[1]],
                 tolerance = 1e-10)
    expect_equal(pairwise$table$p_adjusted, pairwise$table$p_perm)
    expect_equal(maximum$table$p_perm, maximum$table$p_maxT)
    expect_match(contrast$design_type, omnibus$assignment_level)
  }
})

test_that("lower-level validators reject wrong level and confounded actions", {
  child <- .lower_child_example()
  child_design <- .lower_design(child)
  expect_error(
    permanova_unbalanced_nested_lower(
      child_design, level = "observation", permutations = 0
    ),
    "does not match detected"
  )

  balanced <- expand.grid(
    treatment = factor(c("A", "B")), site_local = 1:2, region = factor(1:2)
  )
  balanced$site <- interaction(balanced$region, balanced$site_local)
  balanced$y <- seq_len(nrow(balanced))
  expect_error(
    permanova_unbalanced_nested_lower(
      .lower_design(balanced), permutations = 0
    ),
    "requires an unbalanced hierarchy"
  )

  top <- child
  top$treatment <- factor(ifelse(top$region == levels(top$region)[[1]],
                                 "A", "B"))
  expect_error(
    permanova_unbalanced_nested_lower(.lower_design(top), permutations = 0),
    "top hierarchy level"
  )

  rows <- list()
  site <- 0L
  for (region in 1:2) {
    for (treatment in c("A", "B")) {
      site <- site + 1L
      size <- if (treatment == "A") 2L else 3L
      rows[[length(rows) + 1L]] <- data.frame(
        region = region, site = site, treatment = treatment,
        observation = seq_len(size)
      )
    }
  }
  confounded <- do.call(rbind, rows)
  confounded$region <- factor(confounded$region)
  confounded$site <- factor(confounded$site)
  confounded$treatment <- factor(confounded$treatment)
  confounded$y <- seq_len(nrow(confounded))
  expect_error(
    permanova_unbalanced_nested_lower(
      .lower_design(confounded), permutations = 0
    ),
    "do not span every treatment contrast"
  )

  disconnected <- data.frame(
    region = factor(c(1, 1, 1, 1, 2, 2)),
    site = factor(c(1, 1, 2, 2, 3, 3)),
    treatment = factor(c("A", "B", "A", "B", "C", "C")),
    y = c(-1, 0, -.5, .3, 1, 2)
  )
  expect_error(
    permanova_unbalanced_nested_lower(
      .lower_design(disconnected), permutations = 0
    ),
    "incidence is disconnected"
  )

  no_df <- data.frame(
    region = factor(c(1, 1, 1, 1, 2, 2)),
    site = factor(rep(1:3, each = 2)),
    treatment = factor(rep(c("A", "B", "A"), each = 2)),
    y = c(-1, -.8, .2, .4, 1, 1.2)
  )
  expect_error(
    permanova_unbalanced_nested_lower(.lower_design(no_df), permutations = 0),
    "No experimental-unit residual degrees"
  )

  expect_error(
    permanova_unbalanced_nested_lower(
      .lower_design(.lower_child_example(zero_residual = TRUE)),
      permutations = 0
    ),
    "mean square is numerically zero"
  )
})

test_that("lower-level inference retains semimetric geometry", {
  data <- .lower_observation_example()
  distance <- as.matrix(dist(data$y))
  distance[1, nrow(data)] <- distance[nrow(data), 1] <- 5 * max(distance)
  design <- compile_permanova_design(
    ~ treatment + (1 | region / site), data, distance
  )
  fit <- permanova_unbalanced_nested_lower(
    design, permutations = 15, seed = 4
  )

  expect_false(design$spectral$euclidean)
  expect_true(is.finite(fit$table$p_perm[[1]]))
  expect_match(paste(fit$test$warnings, collapse = " "), "non-Euclidean")
})
