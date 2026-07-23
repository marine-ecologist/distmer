.nested_covariate_example <- function(seed = 7) {
  set.seed(seed)
  d <- expand.grid(rep = 1:5, Area = factor(1:2), Site = factor(1:8))
  d$Location <- factor(c("A", "B", "C", "D")[ceiling(as.integer(d$Site) / 2)])
  site_effect <- rnorm(8, sd = 1.5)[d$Site]
  area_effect <- rnorm(16, sd = 1)[interaction(d$Site, d$Area, drop = TRUE)]
  d$Volume <- rnorm(nrow(d), 80, 20)
  d$y <- 2 * (d$Location == "B") + 0.03 * d$Volume +
    site_effect + area_effect + rnorm(nrow(d))
  d
}

test_that("covariate-adjusted nested numerator equals sequential lm SS", {
  d <- .nested_covariate_example()
  fit <- permanova(
    dist(d$y) ~ Location + Volume + (1 | Site/Area), d, permutations = 0
  )
  table <- fit$fit$table
  location_ref <- anova(lm(y ~ Volume + Location, d))["Location", "Sum Sq"]
  expect_equal(table$ss[table$term == "Location"], location_ref, tolerance = 1e-9)
})

test_that("covariate-adjusted nested Location uses the Site(Location) denominator", {
  d <- .nested_covariate_example()
  base <- permanova_nested(
    compile_permanova_design(~ Location + (1 | Site/Area), d, dist(d$y)),
    permutations = 0
  )
  fit <- permanova_nested_covariate(
    compile_permanova_design(~ Location + Volume + (1 | Site/Area), d, dist(d$y)),
    permutations = 0
  )
  site_ms <- base$table$ms[base$table$term == "Site within fixed model"]
  location_denom_ms <- fit$tests$Location$denominator$ms
  expect_equal(location_denom_ms, site_ms, tolerance = 1e-9)
})

test_that("degrees of freedom partition every observation", {
  d <- .nested_covariate_example()
  fit <- permanova_nested_covariate(
    compile_permanova_design(~ Location + Volume + (1 | Site/Area), d, dist(d$y)),
    permutations = 0
  )
  expect_equal(sum(fit$table$df), nrow(d) - 1L)
})

test_that("variance components are reported over the covariate-adjusted residual", {
  d <- .nested_covariate_example()
  fit <- permanova_nested_covariate(
    compile_permanova_design(~ Location + Volume + (1 | Site/Area), d, dist(d$y)),
    permutations = 0
  )
  components <- fit$variance_components$unrestricted
  expect_named(components, c("residual", "Site", "Site:Area"))
  # the residual component equals the covariate-adjusted within-unit mean square
  within_ms <- fit$table$ms[fit$table$term == "Residual"]
  expect_equal(unname(components[["residual"]]), within_ms, tolerance = 1e-9)
})

test_that("covariate-by-factor interaction (homogeneity of slopes) is tested, not dropped", {
  d <- .nested_covariate_example()
  fit <- permanova_nested_covariate(
    compile_permanova_design(~ Location * Volume + (1 | Site/Area), d, dist(d$y)),
    permutations = 0
  )
  expect_true("Location:Volume" %in% fit$table$term)
  expect_match(fit$tests[["Location:Volume"]]$denominator$expression, "Residual")
  expect_equal(sum(fit$table$df), nrow(d) - 1L)
  # the interaction takes its 3 df from the residual (60 vs 63 without it)
  additive <- permanova_nested_covariate(
    compile_permanova_design(~ Location + Volume + (1 | Site/Area), d, dist(d$y)),
    permutations = 0
  )
  interaction_res <- fit$table$df[fit$table$term == "Residual"]
  additive_res <- additive$table$df[additive$table$term == "Residual"]
  expect_equal(additive_res - interaction_res, 3L)
})

test_that("fixed-effect components are reported separately from random variances", {
  d <- .nested_covariate_example()
  fit <- permanova(dist(d$y) ~ Location + Volume + (1 | Site/Area), d, permutations = 0)
  fixed <- distmer:::.nested_covariate_fixed_components(fit)
  expect_true(all(c("S(Location)", "S(Volume)") %in% fixed$Source))
  # the random variance components stay factor/covariate-free
  expect_named(fit$fit$variance_components$unrestricted,
               c("residual", "Site", "Site:Area"))
})

test_that("the engine requires both a factor and a covariate", {
  d <- .nested_covariate_example()
  factor_only <- compile_permanova_design(
    ~ Location + (1 | Site/Area), d, dist(d$y)
  )
  expect_error(permanova_nested_covariate(factor_only), "continuous covariate")
})

test_that("decompose_covariate splits into zero-mean within and unit-mean between", {
  d <- .nested_covariate_example()
  split <- decompose_covariate(d, "Volume", by = "Site")
  expect_true(all(c("Volume_between", "Volume_within") %in% names(split)))
  # within has mean zero within every Site
  expect_lt(max(abs(tapply(split$Volume_within, split$Site, mean))), 1e-10)
  # between is constant within each Site
  expect_lt(max(tapply(split$Volume_between, split$Site,
                       function(z) diff(range(z)))), 1e-10)
})

test_that("decompose_covariate(center = FALSE) reconstructs the covariate", {
  d <- .nested_covariate_example()
  raw <- decompose_covariate(d, "Volume", by = "Site", center = FALSE)
  expect_equal(raw$Volume_between + raw$Volume_within, d$Volume, tolerance = 1e-12)
})

test_that("decompose_covariate validates its inputs", {
  d <- .nested_covariate_example()
  expect_error(decompose_covariate(d, "nope", by = "Site"), "numeric columns")
  expect_error(decompose_covariate(d, "Volume", by = "nope"), "grouping columns")
  expect_error(decompose_covariate(d, "Location", by = "Site"), "ordinary numeric")
})

test_that("nested covariate Type III adjusts a main slope for its interaction", {
  set.seed(904)
  d <- expand.grid(rep = 1:4, AreaLocal = 1:2, SiteLocal = 1:3,
                   Loc = factor(c("A", "B", "C")))
  d$Site <- interaction(d$Loc, d$SiteLocal, drop = TRUE)
  d$Area <- interaction(d$Site, d$AreaLocal, drop = TRUE)
  loc <- as.integer(d$Loc)
  d$x <- rnorm(nrow(d), mean = c(-1, 0.5, 2)[loc],
               sd = c(0.4, 1, 2)[loc])
  d$y <- 0.7 * d$x + c(-1.2, 0.2, 1.1)[loc] * d$x + rnorm(nrow(d))

  sequential <- permanova(
    dist(d$y) ~ Loc * x + (1 | Site/Area), d,
    type = "I", permutations = 0, mc_permutations = 0
  )
  partial <- permanova(
    dist(d$y) ~ Loc * x + (1 | Site/Area), d,
    type = "III", permutations = 0, mc_permutations = 0
  )
  expect_false(isTRUE(all.equal(
    sequential$fit$tests$x$statistic, partial$fit$tests$x$statistic
  )))

  design <- partial$design
  X <- design$X
  assignment <- attr(X, "assign")
  labels <- attr(terms(design$fixed_formula), "term.labels")
  factor_design <- compile_permanova_design(
    ~ Loc + (1 | Site/Area), design$data, design$distance
  )
  deepest <- distmer:::.validate_nested_design(factor_design)$strata[[3L]]$H
  x_columns <- which(assignment == match("x", labels))
  interaction_columns <- which(assignment == match("Loc:x", labels))
  reduced <- projection_matrix(deepest %*% X[, interaction_columns, drop = FALSE])
  full <- projection_matrix(cbind(
    deepest %*% X[, interaction_columns, drop = FALSE],
    deepest %*% X[, x_columns, drop = FALSE]
  ))
  reference_H <- full$H - reduced$H
  expect_equal(partial$fit$tests$x$numerator$projector, reference_H,
               tolerance = 1e-8)
})

test_that("nested covariate permutation reduced model retains adjusted slopes", {
  d <- .nested_covariate_example()
  fit <- permanova_nested_covariate(
    compile_permanova_design(
      ~ Location * Volume + (1 | Site/Area), d, dist(d$y)
    ), permutations = 0, hypothesis = "I"
  )
  interaction <- fit$tests[["Location:Volume"]]
  main <- fit$tests[["Volume"]]
  expect_equal(
    interaction$reduced_model$projector %*% main$numerator$projector,
    main$numerator$projector, tolerance = 1e-8
  )
})

test_that("nested covariate rejects a structurally higher-level numeric term", {
  d <- .nested_covariate_example()
  d <- decompose_covariate(d, "Volume", by = c("Site", "Area"))
  expect_error(
    permanova(
      dist(d$y) ~ Location + Volume_between + Volume_within +
        (1 | Site/Area), d, permutations = 0, mc_permutations = 0
    ),
    "no estimable variation in the deepest nested stratum"
  )
})
