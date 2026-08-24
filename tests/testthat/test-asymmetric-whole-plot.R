# Whole-plot permutation in the asymmetric engine.
#
# Intact-unit permutation swaps blocks of observation indices, which is only
# defined between units holding the same counts in the same sub-plot classes.
# Real surveys give every unit a distinct count vector, so no two units are
# exchangeable and every whole-plot term returns p = 1 whatever the effect.
# `whole_plot` decides what happens then; these tests pin both branches and,
# critically, that an equally replicated design is untouched.

.asym_fixture <- function(equal, seed = 9L) {
  set.seed(seed)
  cells <- list(c("Hawaii", "P"), c("Hawaii", "U"), c("Mariana", "P"),
                c("Mariana", "U"), c("Samoa", "P"), c("Samoa", "U"),
                c("Line", "U"), c("Marshall", "U"), c("PTNC", "U"))
  rows <- list(); unit <- 0L
  for (cell in cells) {
    for (i in 1:4) {
      unit <- unit + 1L
      for (depth in c("Shallow", "Mid", "Deep")) {
        n <- if (equal) 6L else sample(8:25, 1L)
        for (r in seq_len(n)) {
          rows[[length(rows) + 1L]] <- data.frame(
            ISLAND = sprintf("I%02d", unit), ECOREGION = cell[[1L]],
            POP_STATUS = cell[[2L]], DEPTH = depth, stringsAsFactors = FALSE)
        }
      }
    }
  }
  d <- do.call(rbind, rows)
  for (v in names(d)) d[[v]] <- factor(d[[v]])
  d$DEPTH <- factor(d$DEPTH, levels = c("Shallow", "Mid", "Deep"))
  d$y <- 0.6 * (d$ECOREGION == "Samoa") + 0.4 * (d$DEPTH == "Deep") +
    stats::rnorm(nlevels(d$ISLAND), 0, 0.7)[d$ISLAND] +
    stats::rnorm(nrow(d), 0, 0.6)
  d
}

.asym_design <- function(d) {
  compile_permanova_design(
    ~ DEPTH * ECOREGION * POP_STATUS + (1 | ISLAND) + (1 | ISLAND:DEPTH),
    d, dist(d$y)
  )
}

test_that("an equally replicated design is unchanged by the new argument", {
  design <- .asym_design(.asym_fixture(equal = TRUE))
  auto <- permanova_asymmetric_mixed(design, permutations = 199L, seed = 1L)
  forced <- permanova_asymmetric_mixed(design, permutations = 199L, seed = 1L,
                                       whole_plot = "observation")
  expect_equal(auto$whole_plot, "observation")
  expect_true(auto$whole_plot_exchangeable)
  # bit-identical: nothing already published moves
  expect_equal(auto$table$pseudo_f, forced$table$pseudo_f)
  expect_equal(auto$table$p_perm, forced$table$p_perm)
  # and no estimand note, because no fallback happened
  expect_length(auto$tests$ECOREGION$warnings, 0L)
})

test_that("unequal units make every whole-plot term untestable at observation level", {
  design <- .asym_design(.asym_fixture(equal = FALSE))
  forced <- permanova_asymmetric_mixed(design, permutations = 199L, seed = 1L,
                                       whole_plot = "observation")
  for (term in c("ECOREGION", "POP_STATUS", "ECOREGION:POP_STATUS")) {
    spec <- forced$tests[[term]]
    expect_equal(spec$p_value, 1)                    # only the identity is attainable
    expect_equal(spec$distinct_statistics, 1L)
    expect_equal(spec$minimum_attainable_p, 1)
    expect_true(any(grepl("distinct statistic", spec$warnings)))
  }
})

test_that("the centroid fallback yields real whole-plot p-values and says so", {
  design <- .asym_design(.asym_fixture(equal = FALSE))
  auto <- permanova_asymmetric_mixed(design, permutations = 199L, seed = 1L)
  expect_equal(auto$whole_plot, "centroid")
  expect_false(auto$whole_plot_exchangeable)

  for (term in c("ECOREGION", "POP_STATUS", "ECOREGION:POP_STATUS")) {
    spec <- auto$tests[[term]]
    expect_lt(spec$p_value, 1)
    expect_gt(spec$p_value, 0)
    expect_equal(spec$whole_plot_action, "centroid")
    expect_true(any(grepl("equally weighted", spec$warnings)))
  }
  # sub-plot terms keep the observation-level action and the DEPTH x unit stratum
  depth <- auto$tests$DEPTH
  expect_true(is.na(depth$whole_plot_action))
  expect_match(depth$denominator$expression, "IS(ISLAND)", fixed = TRUE)
})

test_that("the two whole-plot actions target different estimands under imbalance", {
  design <- .asym_design(.asym_fixture(equal = FALSE))
  observation <- permanova_asymmetric_mixed(design, permutations = 0L,
                                            whole_plot = "observation")
  centroid <- permanova_asymmetric_mixed(design, permutations = 0L,
                                         whole_plot = "centroid")
  f_obs <- observation$table$pseudo_f[observation$table$term == "ECOREGION"]
  f_cen <- centroid$table$pseudo_f[centroid$table$term == "ECOREGION"]
  # weighting units by effort is not the same as weighting them equally
  expect_false(isTRUE(all.equal(f_obs, f_cen)))
  # but sub-plot terms are identical either way
  expect_equal(observation$table$pseudo_f[observation$table$term == "DEPTH"],
               centroid$table$pseudo_f[centroid$table$term == "DEPTH"])
})

test_that("forcing centroids on a balanced design reproduces the observation test", {
  design <- .asym_design(.asym_fixture(equal = TRUE))
  observation <- permanova_asymmetric_mixed(design, permutations = 0L,
                                            whole_plot = "observation")
  centroid <- permanova_asymmetric_mixed(design, permutations = 0L,
                                         whole_plot = "centroid")
  for (term in c("ECOREGION", "POP_STATUS", "ECOREGION:POP_STATUS")) {
    expect_equal(observation$table$pseudo_f[observation$table$term == term],
                 centroid$table$pseudo_f[centroid$table$term == term],
                 tolerance = 1e-8, info = term)
  }
})
