.between_plot_example <- function() {
  dat <- expand.grid(replicate = seq_len(3), site = factor(seq_len(6)))
  dat$treatment <- factor(ifelse(as.integer(dat$site) <= 3, "control", "impact"))
  site_effect <- c(-1.2, 0.1, 0.8, -0.5, 0.4, 1.1)
  dat$y <- 2 * (dat$treatment == "impact") +
    site_effect[as.integer(dat$site)] + c(-0.1, 0, 0.1)
  dat
}

test_that("EMS selects the whole-plot residual denominator", {
  dat <- .between_plot_example()
  design <- compile_permanova_design(
    ~ treatment + (1 | site), dat, distance = dist(dat$y)
  )
  ems <- derive_ems(design)

  expect_equal(ems$treatment$numerator_under_null, c(residual = 1, random_intercept = 3))
  expect_equal(ems$treatment$denominator, ems$treatment$numerator_under_null)
  expect_true(ems$treatment$exact_denominator)
  expect_equal(ems$treatment$denominator_df, 4)
})

test_that("Euclidean mixed pseudo-F equals ANOVA on whole-plot means", {
  dat <- .between_plot_example()
  design <- compile_permanova_design(
    ~ treatment + (1 | site), dat, distance = dist(dat$y)
  )
  result <- suppressWarnings(permanova_mixed(
    design, permutations = 0, keep_permutations = TRUE
  ))
  site_means <- aggregate(y ~ site + treatment, dat, mean)
  reference <- summary(aov(y ~ treatment, site_means))[[1]][["F value"]][1]

  expect_equal(result$table$pseudo_f[[1]], unname(reference), tolerance = 1e-10)
  expect_s3_class(result$tests$treatment, "permanova_test_spec")
  expect_match(result$tests$treatment$permutation_algorithm, "Freedman-Lane")
  expect_true(is.na(result$tests$treatment$p_value))
  expect_named(result$tests$treatment, c(
    "term", "hypothesis", "numerator", "denominator", "symbolic_ems",
    "quasi_f_coefficients", "statistic", "p_value", "exchangeable_units",
    "permutation_algorithm", "attainable_permutations", "permutations_used",
    "unique_permutations_used", "exact_permutation", "approximation_status",
    "warnings", "permutation_statistics"
  ))
})

test_that("small attainable spaces are enumerated exactly", {
  dat <- .between_plot_example()
  design <- compile_permanova_design(
    ~ treatment + (1 | site), dat, distance = dist(dat$y)
  )
  result <- suppressWarnings(permanova_mixed(
    design, permutations = 719, seed = 2, keep_permutations = TRUE
  ))

  expect_true(result$permutation$exact)
  expect_equal(result$permutation$attainable, factorial(6))
  expect_equal(result$permutation$used, factorial(6))
  expect_equal(result$tests$treatment$p_value * factorial(6) %% 1, 0, tolerance = 1e-12)
})

test_that("unsupported Phase 2 designs fail explicitly", {
  dat <- .between_plot_example()
  dat$within <- factor(rep(c("a", "b", "c"), times = 6))
  varying <- compile_permanova_design(
    ~ within + (1 | site), dat, distance = dist(dat$y)
  )
  expect_error(permanova_mixed(varying, permutations = 0), "constant within")

  nested <- compile_permanova_design(
    ~ treatment + (1 | treatment/site), dat, distance = dist(dat$y)
  )
  expect_error(permanova_mixed(nested, permutations = 0), "exactly one")

  singleton <- dat[!duplicated(dat$site), ]
  no_replication <- compile_permanova_design(
    ~ treatment + (1 | site), singleton, distance = dist(singleton$y)
  )
  expect_error(permanova_mixed(no_replication, permutations = 0), "At least two")
})

test_that("balanced random-intercept terms honour Type I and Type III hypotheses", {
  skip_if_not_installed("car")
  set.seed(55)
  cells <- data.frame(
    A = factor(c(rep("a1", 7), rep("a2", 7))),
    B = factor(c(rep("b1", 3), rep("b2", 2), rep("b3", 2),
                 rep("b1", 2), rep("b2", 3), rep("b3", 2)))
  )
  cells$Site <- factor(seq_len(nrow(cells)))
  d <- cells[rep(seq_len(nrow(cells)), each = 4), ]
  d$y <- c(0, 1)[d$A] + c(0, 0.4, 1.2)[d$B] +
    c(0, 0.8, -0.3, 0.2, 0, -0.6)[interaction(d$A, d$B)] +
    rnorm(14)[d$Site] + rnorm(nrow(d), sd = 0.3)

  type1 <- permanova(
    dist(d$y) ~ A * B + (1 | Site), d, type = "I",
    permutations = 0, mc_permutations = 0
  )
  type3 <- permanova(
    dist(d$y) ~ A * B + (1 | Site), d, type = "III",
    permutations = 0, mc_permutations = 0
  )
  site <- aggregate(y ~ A + B + Site, d, mean)
  reference1 <- anova(lm(y ~ A * B, site))
  reference3 <- car::Anova(
    lm(y ~ A * B, site, contrasts = list(A = contr.sum, B = contr.sum)),
    type = 3
  )
  expect_equal(type1$table$`Pseudo-F`[type1$table$Source == "A"],
               reference1["A", "F value"], tolerance = 1e-8)
  expect_equal(type3$table$`Pseudo-F`[type3$table$Source == "A"],
               reference3["A", "F value"], tolerance = 1e-8)
  expect_false(isTRUE(all.equal(
    type1$table$`Pseudo-F`[type1$table$Source == "A"],
    type3$table$`Pseudo-F`[type3$table$Source == "A"]
  )))
})
