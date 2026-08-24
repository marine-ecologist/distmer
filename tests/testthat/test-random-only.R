# Models with no fixed term. distmer tests fixed effects against random strata,
# so `permanova()` refuses them -- but the variance-component estimators do
# decompose a pure random structure, and the refusal points there.

test_that("permanova() refuses a model with no fixed term and names the estimator", {
  set.seed(1)
  d <- expand.grid(rep = 1:3, b = factor(1:2), a = factor(1:4))
  d$y <- stats::rnorm(nrow(d))
  D <- dist(d$y)

  expect_error(permanova(~ (1 | a), d, distance = D, permutations = 9L),
               "No fixed term in the model")
  # a `/` chain is pointed at the nested estimator
  expect_error(permanova(~ (1 | a / b), d, distance = D, permutations = 9L),
               "estimate_nested_variance_components")
  # anything else at the general one
  expect_error(permanova(~ (1 | a) + (1 | b), d, distance = D, permutations = 9L),
               "estimate_variance_components")
  # and the message suggests a concrete repair using the user's own grouping
  expect_error(permanova(~ (1 | a), d, distance = D, permutations = 9L),
               "(1 | a)", fixed = TRUE)

  # a model that does have a fixed term is unaffected
  d$tr <- factor(ifelse(as.integer(d$a) <= 2L, "control", "impact"))
  expect_no_error(
    suppressWarnings(permanova(~ tr + (1 | a), d, distance = D, permutations = 9L))
  )
})

test_that("nested variance components do not require a fixed term", {
  set.seed(1)
  d <- expand.grid(rep = 1:3, b = factor(1:2), a = factor(1:4))
  d$y <- stats::rnorm(4)[d$a] +
    0.5 * stats::rnorm(8)[interaction(d$a, d$b)] +
    stats::rnorm(nrow(d), 0, 0.7)

  components <- estimate_nested_variance_components(
    compile_permanova_design(~ (1 | a / b), d, dist(d$y))
  )
  strata <- summary(stats::aov(y ~ 1 + Error(a / b), d))
  within <- strata[["Error: Within"]][[1L]][["Mean Sq"]][[1L]]
  child <- strata[["Error: a:b"]][[1L]][["Mean Sq"]][[1L]]
  parent <- strata[["Error: a"]][[1L]][["Mean Sq"]][[1L]]
  reference <- c(residual = within,
                 `a:b` = (child - within) / 3,
                 a = (parent - child) / 6)

  expect_equal(components$unrestricted[names(reference)], reference,
               tolerance = 1e-10)
})

test_that("the testing path still requires a fixed term", {
  set.seed(1)
  d <- expand.grid(rep = 1:3, b = factor(1:2), a = factor(1:4))
  d$y <- stats::rnorm(nrow(d))
  design <- compile_permanova_design(~ (1 | a / b), d, dist(d$y))
  expect_error(permanova_nested(design, permutations = 0),
               "At least one fixed-effect term is required")
})
