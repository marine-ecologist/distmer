.adjusted_example <- function() {
  data <- expand.grid(
    replicate = 1:3,
    treatment = factor(c("A", "B", "C"), levels = c("A", "B", "C")),
    habitat = factor(c("open", "closed"))
  )
  data <- data[-c(2, 11, 16), ]
  rownames(data) <- NULL
  data$gradient <- seq(-1.3, 1.5, length.out = nrow(data)) +
    c(A = -.2, B = .1, C = .3)[data$treatment]
  data$y <- c(A = -0.4, B = 0.5, C = 1.2)[data$treatment] +
    c(open = -.3, closed = .4)[data$habitat] +
    .8 * data$gradient +
    .25 * (data$treatment == "C") * data$gradient +
    sin(seq_len(nrow(data))) / 10
  data
}

.manual_adjusted_rows <- function(design, at) {
  levels <- levels(design$data$treatment)
  habitat <- levels(design$data$habitat)
  do.call(rbind, lapply(levels, function(level) {
    grid <- expand.grid(
      treatment = factor(level, levels = levels),
      habitat = factor(habitat, levels = habitat),
      gradient = at, KEEP.OUT.ATTRS = FALSE
    )
    colMeans(model.matrix(terms(design$fixed_formula), grid))
  }))
}

test_that("adjusted means reproduce an explicit manual reference grid", {
  data <- .adjusted_example()
  design <- compile_permanova_design(
    ~ treatment * gradient + habitat, data, dist(data$y)
  )
  means <- derive_adjusted_means(
    design, "treatment", at = c(gradient = 0.25), weights = "equal"
  )
  manual <- .manual_adjusted_rows(design, 0.25)

  expect_equal(unname(means$coefficient_rows), unname(manual),
               tolerance = 1e-12)
  expect_equal(means$reference_values, c(gradient = 0.25))
  expect_equal(means$nuisance_weights, c(.5, .5))
})

test_that("adjusted contrast F equals scalar GLH", {
  data <- .adjusted_example()
  design <- compile_permanova_design(
    ~ treatment * gradient + habitat, data, dist(data$y)
  )
  contrast <- c(A = 1, B = -1, C = 0)
  fit <- permanova_adjusted_contrasts(
    design, "treatment", contrast, at = c(gradient = 0.25),
    permutations = 0
  )
  rows <- .manual_adjusted_rows(design, 0.25)
  X <- design$X
  beta <- qr.solve(X, data$y)
  L <- as.numeric(contrast %*% rows)
  estimate <- sum(L * beta)
  numerator <- estimate ^ 2 / as.numeric(L %*% solve(crossprod(X), L))
  denominator <- sum((data$y - X %*% beta)^2) /
    (nrow(X) - qr(X)$rank)

  expect_equal(fit$table$pseudo_f, numerator / denominator,
               tolerance = 1e-10)
  expect_equal(fit$tests[[1]]$reference_values, c(gradient = .25))
})

test_that("adjusted permutations equal direct response reconstruction", {
  data <- .adjusted_example()[1:6, ]
  data$treatment <- factor(rep(c("A", "B", "C"), each = 2))
  data$habitat <- factor(rep(c("open", "closed"), 3))
  data$gradient <- seq(-1, 1, length.out = 6)
  data$y <- c(-1, -.2, .1, 1, 1.3, 2.8)
  design <- compile_permanova_design(
    ~ treatment + gradient, data, dist(data$y)
  )
  contrast <- c(A = 1, B = -1, C = 0)
  fit <- permanova_adjusted_contrasts(
    design, "treatment", contrast, at = c(gradient = 0),
    permutations = "all", keep_permutations = TRUE
  )
  means <- derive_adjusted_means(
    design, "treatment", at = c(gradient = 0)
  )
  L <- as.numeric(contrast %*% means$coefficient_rows)
  projector <- distmer:::.linear_contrast_projector(
    design$X, L, design$projections$fixed$H, 1e-10
  )
  orders <- distmer:::.all_orders(6)
  residual <- diag(6) - projector$H_reduced
  denominator_H <- diag(6) - design$projections$fixed$H
  denominator_df <- 6 - qr(design$X)$rank
  direct <- apply(orders, 1, function(order) {
    y <- (projector$H_reduced + residual[order, , drop = FALSE]) %*% data$y
    sum(y * (projector$H %*% y)) /
      (sum(y * (denominator_H %*% y)) / denominator_df)
  })

  expect_equal(fit$tests[[1]]$permutation_statistics, direct,
               tolerance = 1e-10)
  expect_equal(fit$permutation$used, factorial(6))
})

test_that("adjusted pairwise families and proportional weights are explicit", {
  data <- .adjusted_example()
  design <- compile_permanova_design(
    ~ treatment + habitat + gradient, data, dist(data$y)
  )
  means <- derive_adjusted_means(
    design, "treatment", at = c(gradient = mean(data$gradient)),
    weights = "proportional"
  )
  fit <- permanova_adjusted_pairwise(
    design, "treatment", at = c(gradient = mean(data$gradient)),
    weights = "proportional", permutations = 19, seed = 2
  )

  expect_equal(means$nuisance_weights,
               as.numeric(table(data$habitat)) / nrow(data))
  expect_equal(nrow(fit$table), 3)
  expect_equal(fit$table$p_adjusted, p.adjust(fit$table$p_perm, "holm"))
})

test_that("continuous reference contracts reject implicit or invalid values", {
  data <- .adjusted_example()
  design <- compile_permanova_design(
    ~ treatment + gradient, data, dist(data$y)
  )
  expect_error(derive_adjusted_means(design, "treatment"), "`at` must")
  expect_error(
    derive_adjusted_means(design, "treatment", at = c(wrong = 0)),
    "every continuous nuisance"
  )
  expect_error(
    derive_adjusted_means(design, "treatment", at = c(gradient = Inf)),
    "finite"
  )

  mixed <- compile_permanova_design(
    ~ treatment + gradient + (1 | habitat), data, dist(data$y)
  )
  expect_error(
    derive_adjusted_means(mixed, "treatment", at = c(gradient = 0)),
    "fixed-effect models only"
  )
})

test_that("near-aliased covariates reject before projector construction", {
  set.seed(23)
  treatment <- factor(rep(c("A", "B", "C"), each = 10))
  covariate <- as.numeric(treatment) + 1e-9 * rnorm(length(treatment))
  response <- rnorm(length(treatment))
  design <- compile_permanova_design(
    ~ treatment + covariate,
    data.frame(treatment, covariate, response),
    dist(response)
  )

  expect_error(
    permanova_adjusted_contrasts(
      design, "treatment", c(A = 1, B = -1, C = 0),
      at = c(covariate = 2), permutations = 0
    ),
    "numerically full-rank fixed model;.*ill-conditioned"
  )

  set.seed(23)
  covariate <- as.numeric(treatment) + 1e-2 * rnorm(length(treatment))
  stable <- compile_permanova_design(
    ~ treatment + covariate,
    data.frame(treatment, covariate, response),
    dist(response)
  )
  expect_no_error(derive_adjusted_means(
    stable, "treatment", at = c(covariate = 2)
  ))
  stable_means <- derive_adjusted_means(
    stable, "treatment", at = c(covariate = 2)
  )
  expect_gt(stable_means$fixed_model_condition_ratio, sqrt(1e-10))
})

test_that("adjusted empirical inference retains non-Euclidean geometry", {
  data <- .adjusted_example()
  distance <- as.matrix(dist(data$y))
  distance[1, nrow(data)] <- distance[nrow(data), 1] <- 5 * max(distance)
  design <- compile_permanova_design(
    ~ treatment + gradient, data, distance
  )
  fit <- permanova_adjusted_contrasts(
    design, "treatment", c(A = 1, B = -1, C = 0),
    at = c(gradient = 0), permutations = 9, seed = 3
  )

  expect_false(design$spectral$euclidean)
  expect_true(is.finite(fit$table$p_perm))
  expect_match(fit$tests[[1]]$warnings[[1]], "non-Euclidean")
})
