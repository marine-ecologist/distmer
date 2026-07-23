.unbalanced_additive_example <- function() {
  cells <- expand.grid(A = factor(c("a1", "a2")), B = factor(c("b1", "b2", "b3")))
  counts <- c(2L, 4L, 3L, 2L, 5L, 3L)
  dat <- cells[rep(seq_len(nrow(cells)), counts), , drop = FALSE]
  rownames(dat) <- NULL
  set.seed(601)
  dat$y <- 1.2 * (dat$A == "a2") +
    c(b1 = -0.5, b2 = 0.2, b3 = 0.8)[dat$B] +
    rnorm(nrow(dat), sd = 0.4)
  dat
}

.unbalanced_factorial_example <- function() {
  cells <- expand.grid(A = factor(c("a1", "a2")), B = factor(c("b1", "b2")))
  counts <- c(2L, 4L, 3L, 5L)
  dat <- cells[rep(seq_len(nrow(cells)), counts), , drop = FALSE]
  rownames(dat) <- NULL
  contrasts(dat$A) <- contr.sum(2)
  contrasts(dat$B) <- contr.sum(2)
  set.seed(602)
  dat$y <- 0.8 * (dat$A == "a2") + 0.5 * (dat$B == "b2") +
    1.1 * (dat$A == "a2" & dat$B == "b2") +
    rnorm(nrow(dat), sd = 0.35)
  dat
}

.extra_ss_f <- function(full, reduced) {
  full_residual <- residuals(full)
  reduced_residual <- residuals(reduced)
  numerator_df <- df.residual(reduced) - df.residual(full)
  ((sum(reduced_residual ^ 2) - sum(full_residual ^ 2)) / numerator_df) /
    (sum(full_residual ^ 2) / df.residual(full))
}

test_that("unbalanced Type I hypotheses equal sequential linear-model ANOVA", {
  dat <- .unbalanced_additive_example()
  design <- compile_permanova_design(~ A + B, dat, dist(dat$y))
  fit <- permanova_unbalanced(
    design, hypothesis = "I", permutations = 0
  )
  reference <- summary(aov(y ~ A + B, dat))[[1]]

  expect_s3_class(fit, "permanova_unbalanced")
  expect_equal(
    fit$tests$A$statistic,
    unname(reference["A", "F value"]), tolerance = 1e-10
  )
  expect_equal(
    fit$tests$B$statistic,
    unname(reference["B", "F value"]), tolerance = 1e-10
  )
  expect_equal(fit$hypothesis, "I")
  expect_true(all(vapply(
    fit$tests, function(test) test$estimability$fully_estimable, logical(1)
  )))
})

test_that("unbalanced additive Type II tests are order invariant", {
  dat <- .unbalanced_additive_example()
  forward <- permanova_unbalanced(
    compile_permanova_design(~ A + B, dat, dist(dat$y)),
    hypothesis = "II", permutations = 0
  )
  reverse <- permanova_unbalanced(
    compile_permanova_design(~ B + A, dat, dist(dat$y)),
    hypothesis = "II", permutations = 0
  )
  full <- lm(y ~ A + B, dat)
  reference_A <- .extra_ss_f(full, lm(y ~ B, dat))
  reference_B <- .extra_ss_f(full, lm(y ~ A, dat))

  expect_equal(forward$tests$A$statistic, reference_A, tolerance = 1e-10)
  expect_equal(forward$tests$B$statistic, reference_B, tolerance = 1e-10)
  expect_equal(
    forward$tests$A$statistic, reverse$tests$A$statistic,
    tolerance = 1e-10
  )
  expect_equal(
    forward$tests$B$statistic, reverse$tests$B$statistic,
    tolerance = 1e-10
  )
})

test_that("unbalanced Type III tests match full-versus-drop references", {
  dat <- .unbalanced_factorial_example()
  design <- compile_permanova_design(~ A * B, dat, dist(dat$y))
  fit <- permanova_unbalanced(
    design, hypothesis = "III", permutations = 0
  )
  X <- model.matrix(~ A * B, dat)
  assignment <- attr(X, "assign")
  full <- lm.fit(X, dat$y)
  full_rss <- sum(full$residuals ^ 2)
  full_df <- nrow(X) - full$rank
  references <- setNames(vapply(seq_along(design$fixed_effects), function(k) {
    reduced <- lm.fit(X[, assignment != k, drop = FALSE], dat$y)
    numerator_df <- full$rank - reduced$rank
    ((sum(reduced$residuals ^ 2) - full_rss) / numerator_df) /
      (full_rss / full_df)
  }, numeric(1)), design$fixed_effects)

  expect_equal(
    vapply(fit$tests, `[[`, numeric(1), "statistic"),
    references, tolerance = 1e-10
  )
  expect_equal(fit$tests$A$adjusted_for, c("B", "A:B"))
  expect_equal(fit$tests$`A:B`$adjusted_for, c("A", "B"))
})

test_that("Type III interactions require explicit sum contrasts", {
  dat <- .unbalanced_factorial_example()
  contrasts(dat$A) <- contr.treatment(2)
  contrasts(dat$B) <- contr.treatment(2)
  design <- compile_permanova_design(~ A * B, dat, dist(dat$y))

  expect_error(
    permanova_unbalanced(design, hypothesis = "III", permutations = 0),
    "sum-to-zero contrasts such as `contr.sum`"
  )
})

test_that("incomplete additive models retain estimable main effects", {
  dat <- .unbalanced_additive_example()
  incomplete <- subset(dat, !(A == "a2" & B == "b3"))
  incomplete <- droplevels(incomplete)
  design <- compile_permanova_design(
    ~ A + B, incomplete, dist(incomplete$y)
  )
  fit <- permanova_unbalanced(
    design, hypothesis = "II", permutations = 0
  )

  expect_true(all(vapply(
    fit$hypotheses, `[[`, logical(1), "fully_estimable"
  )))
  expect_true(all(is.finite(vapply(
    fit$tests, `[[`, numeric(1), "statistic"
  ))))
})

test_that("aliased and partially estimable hypotheses reject explicitly", {
  dat <- .unbalanced_additive_example()
  incomplete <- subset(dat, !(A == "a2" & B == "b3"))
  incomplete <- droplevels(incomplete)
  contrasts(incomplete$A) <- contr.sum(2)
  contrasts(incomplete$B) <- contr.sum(3)
  design <- compile_permanova_design(
    ~ A * B, incomplete, dist(incomplete$y)
  )

  expect_error(
    permanova_unbalanced(design, hypothesis = "III", permutations = 0),
    "full-rank fixed model"
  )
  hypotheses <- derive_unbalanced_hypotheses(
    design, hypothesis = "I", terms = "A:B"
  )
  expect_false(hypotheses$`A:B`$fully_estimable)
  expect_match(hypotheses$`A:B`$reason, "estimable df")
  expect_error(
    permanova_unbalanced(
      design, hypothesis = "I", terms = "A:B", permutations = 0
    ),
    "not fully estimable"
  )
})

test_that("Type II main effects respect marginality (adjust for B, exclude A:B)", {
  dat <- .unbalanced_factorial_example()
  design <- compile_permanova_design(~ A * B, dat, dist(dat$y))
  hypotheses <- derive_unbalanced_hypotheses(
    design, hypothesis = "II", terms = "A"
  )
  expect_true(hypotheses$A$fully_estimable)
  expect_identical(hypotheses$A$adjusted_for, "B")
  expect_identical(hypotheses$A$excluded_higher_order, "A:B")
  expect_equal(hypotheses$A$df, hypotheses$A$nominal_df)

  fit <- permanova_unbalanced(design, hypothesis = "II", permutations = 0)
  expect_true(all(c("A", "B", "A:B") %in% fit$table$term))
})

test_that("Type II factor SS match car::Anova(type = 2) on a 1-D response", {
  skip_if_not_installed("car")
  set.seed(1)
  dat <- expand.grid(A = factor(c("a1", "a2", "a3")), B = factor(c("b1", "b2")))
  dat <- dat[rep(seq_len(nrow(dat)), c(4, 3, 5, 2, 6, 3)), ]
  rownames(dat) <- NULL
  dat$y <- rnorm(nrow(dat))
  design <- compile_permanova_design(~ A * B, dat, dist(dat$y))
  fit <- permanova_unbalanced(design, hypothesis = "II", permutations = 0)
  reference <- car::Anova(lm(y ~ A * B, dat), type = 2)
  for (term in c("A", "B", "A:B")) {
    expect_equal(
      fit$table$ss[fit$table$term == term],
      unname(reference[["Sum Sq"]][trimws(rownames(reference)) == term]),
      tolerance = 1e-8
    )
  }
})

test_that("exact unbalanced Freedman-Lane statistics match direct ANOVA", {
  dat <- data.frame(
    group = factor(c("a", "a", "b", "b", "b")),
    y = c(-1.1, 0.3, 0.4, 1.8, 3.2)
  )
  design <- compile_permanova_design(~ group, dat, dist(dat$y))
  fit <- suppressWarnings(permanova_unbalanced(
    design, hypothesis = "III", permutations = "all",
    keep_permutations = TRUE, max_exact = 1000L
  ))
  orders <- distmer:::.all_orders(nrow(dat))
  reduced <- lm(y ~ 1, dat)
  fitted_reduced <- fitted(reduced)
  residual_reduced <- residuals(reduced)
  direct <- apply(orders, 1L, function(index) {
    response <- fitted_reduced + residual_reduced[index]
    unname(summary(aov(response ~ dat$group))[[1]][1L, "F value"])
  })

  expect_true(fit$tests$group$exact_permutation)
  expect_equal(fit$tests$group$attainable_permutations, factorial(5))
  expect_equal(
    fit$tests$group$permutation_statistics,
    direct, tolerance = 1e-10
  )
  expect_equal(
    fit$tests$group$p_value * factorial(5) %% 1,
    0, tolerance = 1e-12
  )
})

test_that("unbalanced non-Euclidean inference retains empirical p-values", {
  dat <- .unbalanced_additive_example()
  distance <- as.matrix(dist(cbind(seq_len(nrow(dat)), dat$y)))
  distance[1L, 2L] <- distance[2L, 1L] <- 3 * max(distance)
  design <- compile_permanova_design(~ A + B, dat, distance)
  fit <- permanova_unbalanced(
    design, hypothesis = "II", permutations = 19, seed = 97
  )

  expect_false(design$spectral$euclidean)
  expect_true(all(is.finite(vapply(
    fit$tests, `[[`, numeric(1), "p_value"
  ))))
  expect_match(
    paste(fit$tests$A$warnings, collapse = " "),
    "complete non-Euclidean"
  )
})

test_that("unbalanced mixed designs and invalid term requests fail", {
  dat <- .unbalanced_additive_example()
  dat$cluster <- factor(rep(seq_len(ceiling(nrow(dat) / 2)), each = 2)[
    seq_len(nrow(dat))
  ])
  mixed <- compile_permanova_design(
    ~ A + (1 | cluster), dat, dist(dat$y)
  )
  expect_error(
    permanova_unbalanced(mixed, hypothesis = "I", permutations = 0),
    "fixed effects only"
  )
  fixed <- compile_permanova_design(~ A + B, dat, dist(dat$y))
  expect_error(
    permanova_unbalanced(
      fixed, hypothesis = "II", terms = "missing", permutations = 0
    ),
    "term labels"
  )
})
