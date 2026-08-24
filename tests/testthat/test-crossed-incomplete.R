# Crossed mixed inference where cells may be empty and replication unequal.
#
# The correctness anchor is agreement: on a complete, equally replicated design
# both this engine and the validated `permanova_crossed()` apply, and they must
# return the same numbers. Everything else this engine does is the same
# construction carried into designs the balanced one refuses.

.crossed_example <- function(seed = 11) {
  set.seed(seed)
  dat <- expand.grid(rep = 1:3, block = factor(1:5), a = factor(c("a1", "a2")),
                     b = factor(c("b1", "b2", "b3")))
  dat$y <- rnorm(nrow(dat)) + as.numeric(dat$a) + as.numeric(dat$block) / 2 +
    as.numeric(dat$a) * as.numeric(dat$block) / 3
  dat
}

.crossed_formula <- ~ a * b + (1 | block) + (1 | a:block) + (1 | b:block) +
  (1 | a:b:block)

test_that("on a complete equally replicated design it matches permanova_crossed", {
  dat <- .crossed_example()
  design <- compile_permanova_design(.crossed_formula, dat, dist(dat$y))
  old <- permanova_crossed(design, permutations = 0)
  new <- permanova_crossed_incomplete(design, permutations = 0)

  for (term in old$table[[1L]]) {
    i <- match(term, old$table[[1L]])
    j <- match(term, new$table$term)
    expect_false(is.na(j), info = term)
    expect_equal(new$table$df[[j]], old$table$df[[i]], info = term)
    expect_equal(new$table$ss[[j]], old$table$ss[[i]], tolerance = 1e-10,
                 info = term)
    if (!is.na(old$table$pseudo_f[[i]])) {
      expect_equal(new$table$pseudo_f[[j]], old$table$pseudo_f[[i]],
                   tolerance = 1e-8, info = term)
    }
  }
})

test_that("an incomplete design reports estimable df, not nominal df", {
  dat <- .crossed_example()
  # Drop two whole a:block cells, so the a:block interaction can no longer
  # carry its nominal (2-1)*(5-1) = 4 dimensions.
  gone <- dat$a == "a2" & dat$block %in% c("4", "5")
  dat <- droplevels(dat[!gone, ])
  design <- compile_permanova_design(.crossed_formula, dat, dist(dat$y))
  fit <- permanova_crossed_incomplete(design, permutations = 0)

  row <- fit$table[fit$table$term == "a:block", ]
  expect_equal(row$nominal_df, 4L)
  expect_lt(row$df, row$nominal_df)
  # occupied a:block cells - a levels - block levels + 1
  occupied <- sum(table(dat$a, dat$block) > 0)
  expect_equal(row$df, occupied - nlevels(dat$a) - nlevels(dat$block) + 1L)
  expect_true("a:block" %in% fit$rank_deficient)

  # the balanced engine cannot do this at all
  expect_error(permanova_crossed(design, permutations = 0),
               "complete, equally replicated")
  # df still add up
  model <- fit$table[!fit$table$term %in% c("Residual", "Total"), ]
  expect_equal(sum(model$df) + fit$table$df[fit$table$term == "Residual"],
               fit$table$df[fit$table$term == "Total"])
})

test_that("unequal replication alone is accepted and totals still balance", {
  dat <- .crossed_example()
  dat <- dat[-c(1, 2, 40), ]   # every cell still present, sizes now differ
  design <- compile_permanova_design(.crossed_formula, dat, dist(dat$y))
  fit <- permanova_crossed_incomplete(design, permutations = 0)
  expect_length(fit$rank_deficient, 0L)
  model <- fit$table[!fit$table$term %in% c("Residual", "Total"), ]
  expect_equal(sum(model$df) + fit$table$df[fit$table$term == "Residual"],
               fit$table$df[fit$table$term == "Total"])
  expect_true(all(fit$table$ss >= -1e-8, na.rm = TRUE))
})

test_that("denominators come from the restricted expected mean squares", {
  dat <- .crossed_example()
  design <- compile_permanova_design(.crossed_formula, dat, dist(dat$y))
  fit <- permanova_crossed_incomplete(design, permutations = 0)

  # A random stratum carries no component from a random-by-fixed interaction
  # sitting inside it: E(MS_block) holds nothing of a:block.
  expect_equal(fit$ems["block", "a:block"], 0)
  expect_equal(fit$ems["block", "a:b:block"], 0)
  # but the interaction stratum does carry its own component
  expect_gt(fit$ems["a:block", "a:block"], 0)

  # A fixed term is tested over its own interaction with the random factor,
  # never over the residual, and never over itself.
  expect_match(fit$table$denominator[fit$table$term == "a"], "a:block")
  expect_match(fit$table$denominator[fit$table$term == "a:b"], "a:b:block")
  # A random term is not its own denominator.
  block <- fit$table$denominator[fit$table$term == "block"]
  expect_false(grepl("^[0-9.]+\\*block", block))
})

test_that("permutation produces usable p-values and honours the seed", {
  dat <- .crossed_example()
  design <- compile_permanova_design(.crossed_formula, dat, dist(dat$y))
  first  <- permanova_crossed_incomplete(design, permutations = 99, seed = 3)
  second <- permanova_crossed_incomplete(design, permutations = 99, seed = 3)
  expect_equal(first$table$p_perm, second$table$p_perm)

  tested <- !is.na(first$table$p_perm)
  expect_gt(sum(tested), 0L)
  expect_true(all(first$table$p_perm[tested] > 0 &
                    first$table$p_perm[tested] <= 1))
  # `a` carries a real effect in the simulation
  expect_lt(first$table$p_perm[first$table$term == "a"], 0.2)
})

test_that("it refuses designs it does not cover", {
  dat <- .crossed_example()
  # A model with no random term is no longer refused -- it is the incomplete
  # fixed case, tested over the residual.
  fixed_only <- permanova_crossed_incomplete(
    compile_permanova_design(~ a * b, dat, dist(dat$y)), permutations = 0)
  expect_length(fixed_only$random_terms, 0L)
  expect_error(
    permanova_crossed_incomplete(
      compile_permanova_design(~ a + (1 | block/rep), dat, dist(dat$y)),
      permutations = 0),
    "nesting")
  expect_error(permanova_crossed_incomplete(dat), "compile_permanova_design")
})

test_that("on a balanced complete design all three hypotheses coincide", {
  dat <- .crossed_example()
  design <- compile_permanova_design(.crossed_formula, dat, dist(dat$y))
  fits <- lapply(c("I", "II", "III"), function(h) {
    permanova_crossed_incomplete(design, hypothesis = h, permutations = 0)
  })
  names(fits) <- c("I", "II", "III")
  # Orthogonality makes the adjustment set irrelevant when every cell is
  # present and equally replicated.
  for (h in c("I", "III")) {
    expect_equal(fits[[h]]$table$df, fits$II$table$df, info = h)
    expect_equal(fits[[h]]$table$ss, fits$II$table$ss, tolerance = 1e-8,
                 info = h)
  }
})

test_that("Type III can lose estimable dimensions that Type II keeps", {
  dat <- .crossed_example()
  gone <- dat$a == "a2" & dat$block %in% c("4", "5")
  dat <- droplevels(dat[!gone, ])
  design <- compile_permanova_design(.crossed_formula, dat, dist(dat$y))
  two   <- permanova_crossed_incomplete(design, hypothesis = "II", permutations = 0)
  three <- permanova_crossed_incomplete(design, hypothesis = "III", permutations = 0)

  # With a term containing it still in the model, a Type III main effect can be
  # spanned by that term and report nothing. Type II compares against the terms
  # that do NOT contain it, so it stays well posed. This is why Type II is the
  # default, and why Type III is offered but not relied on.
  expect_true(all(three$table$df <= two$table$df))
  expect_lt(sum(three$table$df[three$table$term %in% c("a", "b", "a:b")]),
            sum(two$table$df[two$table$term %in% c("a", "b", "a:b")]))
  expect_gt(two$table$df[two$table$term == "a"], 0L)

  expect_error(permanova_crossed_incomplete(design, hypothesis = "IV"),
               "should be one of")
})

test_that("an incomplete fixed design needs no random term", {
  # Two crossed fixed factors where one combination was never sampled. The
  # Type III fixed engine refuses this ("aliased columns from incomplete
  # cells"); the rank-based construction here reports what is estimable.
  set.seed(11)
  dat <- expand.grid(a = factor(c("a1", "a2", "a3")),
                     b = factor(c("b1", "b2", "b3")),
                     rep = 1:4)
  # one cell never sampled: eight of nine remain, so the interaction spans
  # 8 - 3 - 3 + 1 = 3 dimensions of the 4 a complete design would give
  dat <- dat[!(dat$a == "a3" & dat$b == "b3"), ]
  dat$y <- rnorm(nrow(dat)) + as.numeric(dat$a)
  design <- compile_permanova_design(~ a * b, dat, dist(dat$y))

  expect_error(permanova(dist(dat$y) ~ a * b, dat, permutations = 0,
                         type = "III"),
               "full-rank|aliased")

  fit <- permanova_crossed_incomplete(design, permutations = 0)
  expect_s3_class(fit, "permanova_crossed_incomplete")
  expect_length(fit$random_terms, 0L)
  # With no random term every stratum is the residual, so that is the
  # denominator for all of them.
  tested <- fit$table[!is.na(fit$table$denominator), ]
  expect_true(all(tested$denominator == "1*Residual"))
  # a:b has one estimable dimension of the two a complete design would give
  interaction <- fit$table[fit$table$term == "a:b", ]
  expect_equal(interaction$nominal_df, 4L)
  expect_equal(interaction$df, 3L)
  expect_true("a:b" %in% fit$rank_deficient)
  # Type II sums of squares are adjusted, not orthogonal, so they are not a
  # partition of the total; the sequential ones are.
  sequential <- permanova_crossed_incomplete(design, hypothesis = "I",
                                             permutations = 0)
  parts <- sequential$table[sequential$table$term != "Total", ]
  expect_equal(sum(parts$ss),
               sequential$table$ss[sequential$table$term == "Total"],
               tolerance = 1e-8)
})

test_that("a continuous covariate keeps its single degree of freedom", {
  set.seed(21)
  dat <- expand.grid(a = factor(c("a1", "a2")), block = factor(1:5), rep = 1:3)
  dat$depth <- round(runif(nrow(dat), 0, 50), 2)
  dat$y <- 0.4 * dat$depth + as.numeric(dat$a) + rnorm(nrow(dat))
  design <- compile_permanova_design(~ depth + a + (1 | block), dat, dist(dat$y))
  fit <- permanova_crossed_incomplete(design, permutations = 0)

  # A slope is one dimension, not one per distinct value.
  covariate <- fit$table[fit$table$term == "depth", ]
  expect_equal(covariate$df, 1L)
  expect_equal(covariate$nominal_df, 1L)
  expect_equal(fit$table$df[fit$table$term == "a"], 1L)

  # and it explains what a direct regression says it does: with the covariate
  # fitted last, its sum of squares is the drop in residual sum of squares.
  full <- lm(y ~ a + block + depth, dat)
  without <- lm(y ~ a + block, dat)
  expect_equal(fit$table$ss[fit$table$term == "depth"],
               sum(residuals(without)^2) - sum(residuals(full)^2),
               tolerance = 1e-6)
})
