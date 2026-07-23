.contrast_example <- function() {
  data.frame(
    group = factor(c("A", "A", "B", "B", "C"),
                   levels = c("A", "B", "C")),
    y = c(-1.2, 0.2, 0.4, 1.8, 2.5)
  )
}

.direct_contrast_f <- function(y, group, coefficients) {
  means <- tapply(y, group, mean)
  counts <- table(group)
  contrast_ss <- sum(coefficients * means) ^ 2 /
    sum(coefficients ^ 2 / counts)
  full <- lm(y ~ group)
  contrast_ss / (sum(residuals(full) ^ 2) / df.residual(full))
}

test_that("planned contrast F statistics equal scalar contrast formulas", {
  dat <- .contrast_example()
  design <- compile_permanova_design(~ group, dat, dist(dat$y))
  C <- rbind(
    `A - B` = c(A = 1, B = -1, C = 0),
    `AB - C` = c(A = 0.5, B = 0.5, C = -1)
  )
  fit <- suppressWarnings(permanova_contrasts(
    design, C, permutations = "all", keep_permutations = TRUE
  ))
  direct <- apply(C, 1L, function(coefficients) {
    .direct_contrast_f(dat$y, dat$group, coefficients)
  })
  expect_equal(fit$table$pseudo_f, unname(direct), tolerance = 1e-10)
  expect_equal(fit$permutation$attainable, factorial(5))
  expect_equal(fit$table$p_adjusted,
               p.adjust(fit$table$p_perm, method = "holm"))
})

test_that("every exact contrast permutation matches a reduced scalar fit", {
  dat <- .contrast_example()
  design <- compile_permanova_design(~ group, dat, dist(dat$y))
  coefficients <- c(A = 1, B = -1, C = 0)
  fit <- suppressWarnings(permanova_contrasts(
    design, coefficients, permutations = "all", keep_permutations = TRUE
  ))
  X_group <- model.matrix(~ 0 + group, dat)
  Q <- qr.Q(qr(matrix(coefficients, ncol = 1L)), complete = TRUE)
  X_reduced <- X_group %*% Q[, -1L, drop = FALSE]
  reduced <- lm.fit(X_reduced, dat$y)
  fitted_reduced <- as.vector(X_reduced %*% reduced$coefficients)
  residual_reduced <- dat$y - fitted_reduced
  orders <- .all_orders(nrow(dat))
  direct <- apply(orders, 1L, function(index) {
    response <- fitted_reduced + residual_reduced[index]
    .direct_contrast_f(response, dat$group, coefficients)
  })
  expect_equal(fit$tests[[1L]]$permutation_statistics, direct,
               tolerance = 1e-10)
})

test_that("pairwise helper constructs and adjusts the complete family", {
  dat <- .contrast_example()
  design <- compile_permanova_design(~ group, dat, dist(dat$y))
  fit <- suppressWarnings(permanova_pairwise(
    design, permutations = "all", adjust = "bonferroni"
  ))
  expect_equal(fit$table$contrast, c("A - B", "A - C", "B - C"))
  expect_equal(fit$table$p_adjusted,
               p.adjust(fit$table$p_perm, method = "bonferroni"))
  expect_equal(unname(rowSums(fit$contrast_matrix)), c(0, 0, 0))
  expect_match(fit$family, "all pairwise")
})

test_that("Phase 13 rejects unsupported contrast definitions", {
  dat <- .contrast_example()
  design <- compile_permanova_design(~ group, dat, dist(dat$y))
  expect_error(
    permanova_contrasts(design, c(A = 1, B = 0, C = 0), permutations = 0),
    "sum to zero"
  )
  expect_error(
    permanova_contrasts(design, c(A = 1, B = -1), permutations = 0),
    "one coefficient"
  )

  dat$site <- factor(seq_len(nrow(dat)))
  mixed <- compile_permanova_design(
    ~ group + (1 | site), dat, dist(dat$y)
  )
  expect_error(permanova_pairwise(mixed, permutations = 0),
               "fixed one-factor designs only")

  dat$x <- seq_len(nrow(dat))
  multifactor <- compile_permanova_design(~ group + x, dat, dist(dat$y))
  expect_error(permanova_pairwise(multifactor, permutations = 0),
               "exactly one categorical")
})

test_that("Phase 13 retains non-Euclidean empirical inference", {
  dat <- .contrast_example()
  D <- as.matrix(dist(cbind(dat$y, seq_len(nrow(dat)))))
  D[1L, 2L] <- D[2L, 1L] <- 5 * max(D)
  design <- compile_permanova_design(~ group, dat, D)
  fit <- suppressWarnings(permanova_pairwise(
    design, permutations = "all"
  ))
  expect_false(design$spectral$euclidean)
  expect_true(all(is.finite(fit$table$p_perm)))
  expect_true(all(vapply(fit$tests, function(x) {
    any(grepl("non-Euclidean", x$warnings))
  }, logical(1))))
})
