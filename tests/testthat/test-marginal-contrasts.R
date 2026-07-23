.marginal_example <- function() {
  cells <- expand.grid(
    A = factor(c("a1", "a2", "a3")),
    B = factor(c("b1", "b2"))
  )
  counts <- c(2L, 4L, 3L, 5L, 2L, 4L)
  dat <- cells[rep(seq_len(nrow(cells)), counts), ]
  rownames(dat) <- NULL
  dat$y <- c(a1 = -0.4, a2 = 0.5, a3 = 1.1)[dat$A] +
    c(b1 = -0.3, b2 = 0.6)[dat$B] +
    c(-0.2, 0.1, 0.4, 0.3, -0.1, 0.2)[
      interaction(dat$A, dat$B)
    ] + sin(seq_len(nrow(dat))) / 3
  dat
}

.scalar_linear_f <- function(model, coefficient_row) {
  beta <- coef(model)
  information_inverse <- solve(crossprod(model.matrix(model)))
  ss <- as.numeric(coefficient_row %*% beta) ^ 2 /
    as.numeric(coefficient_row %*% information_inverse %*% coefficient_row)
  ss / (sum(residuals(model) ^ 2) / df.residual(model))
}

test_that("equal and proportional marginal means match manual predictions", {
  dat <- .marginal_example()
  design <- compile_permanova_design(~ A * B, dat, dist(dat$y))
  equal <- derive_marginal_means(design, "A", "equal")
  proportional <- derive_marginal_means(design, "A", "proportional")
  model <- lm(y ~ A * B, dat)

  equal_grid <- expand.grid(A = levels(dat$A), B = levels(dat$B))
  equal_reference <- tapply(
    predict(model, equal_grid), equal_grid$A, mean
  )
  B_weights <- prop.table(table(dat$B))
  proportional_reference <- vapply(levels(dat$A), function(level) {
    grid <- data.frame(
      A = factor(level, levels = levels(dat$A)),
      B = factor(levels(dat$B), levels = levels(dat$B))
    )
    sum(predict(model, grid) * B_weights)
  }, numeric(1))
  expect_equal(
    as.vector(equal$coefficient_rows %*% coef(model)),
    as.numeric(equal_reference), tolerance = 1e-10
  )
  expect_equal(
    as.vector(proportional$coefficient_rows %*% coef(model)),
    as.numeric(proportional_reference), tolerance = 1e-10
  )
})

test_that("marginal contrast F and permutations match scalar references", {
  dat <- .marginal_example()
  design <- compile_permanova_design(~ A * B, dat, dist(dat$y))
  contrast <- c(a1 = 1, a2 = -1, a3 = 0)
  fit <- permanova_marginal_contrasts(
    design, "A", contrast, weights = "equal",
    permutations = 19, seed = 42, keep_permutations = TRUE
  )
  model <- lm(y ~ A * B, dat)
  coefficient_row <- fit$tests[[1L]]$marginal_coefficient_row
  direct_f <- .scalar_linear_f(model, coefficient_row)
  expect_equal(fit$tests[[1L]]$statistic, direct_f, tolerance = 1e-10)

  X <- model.matrix(model)
  projector <- .linear_contrast_projector(
    X, coefficient_row, projection_matrix(X)$H, 1e-10
  )
  reduced_y <- projector$H_reduced %*% dat$y
  residual_y <- dat$y - reduced_y
  set.seed(42)
  plan <- .permutation_orders(nrow(dat), 19, 100000L)
  direct <- apply(plan$orders, 1L, function(index) {
    y <- as.vector(reduced_y + residual_y[index])
    scalar <- lm.fit(X, y)
    beta <- scalar$coefficients
    ss <- as.numeric(coefficient_row %*% beta) ^ 2 /
      as.numeric(coefficient_row %*% solve(crossprod(X)) %*% coefficient_row)
    ss / (sum(scalar$residuals ^ 2) / scalar$df.residual)
  })
  expect_equal(fit$tests[[1L]]$permutation_statistics, direct,
               tolerance = 1e-10)
})

test_that("marginal pairwise constructs complete adjusted family", {
  dat <- .marginal_example()
  design <- compile_permanova_design(~ A + B, dat, dist(dat$y))
  fit <- permanova_marginal_pairwise(
    design, "A", weights = "proportional", permutations = 19, seed = 1
  )
  expect_equal(fit$table$contrast, c("a1 - a2", "a1 - a3", "a2 - a3"))
  expect_equal(fit$table$p_adjusted,
               p.adjust(fit$table$p_perm, method = "holm"))
})

test_that("marginal contrasts reject undefined reference grids", {
  dat <- .marginal_example()
  dat$x <- seq_len(nrow(dat))
  continuous <- compile_permanova_design(~ A + x, dat, dist(dat$y))
  expect_error(derive_marginal_means(continuous, "A"),
               "every fixed predictor")

  mixed <- compile_permanova_design(
    ~ A + B + (1 | A), dat, dist(dat$y)
  )
  expect_error(derive_marginal_means(mixed, "A"), "fixed-effect models only")

  incomplete <- dat[!(dat$A == "a3" & dat$B == "b2"), ]
  incomplete <- droplevels(incomplete)
  aliased <- compile_permanova_design(~ A * B, incomplete, dist(incomplete$y))
  expect_error(derive_marginal_means(aliased, "A"), "full-rank")
  expect_error(derive_marginal_means(
    compile_permanova_design(~ A + B, dat, dist(dat$y)), "B:A"
  ), "must name")
})

test_that("marginal contrasts retain non-Euclidean inference", {
  dat <- .marginal_example()
  D <- as.matrix(dist(cbind(dat$y, seq_len(nrow(dat)))))
  D[1L, 2L] <- D[2L, 1L] <- 5 * max(D)
  design <- compile_permanova_design(~ A + B, dat, D)
  fit <- permanova_marginal_pairwise(
    design, "A", permutations = 19, seed = 2
  )
  expect_false(design$spectral$euclidean)
  expect_true(all(is.finite(fit$table$p_perm)))
})
