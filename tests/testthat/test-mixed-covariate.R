# A single-random-intercept ANCOVA is a split-plot: the between-unit factor is
# the whole-plot effect (tested over the group stratum) and the observation-level
# covariate is a sub-plot nuisance (tested over the within-unit residual). The
# correctness anchor is the univariate split-plot model
# aov(y ~ covariate + factor + Error(group)); the covariate is deliberately given
# between-unit structure so the two-slope adjustment is exercised.

.mixed_cov_example <- function(seed = 7) {
  set.seed(seed)
  d <- expand.grid(rep = 1:5, Site = factor(1:12))
  d$Loc <- factor(rep(c("A", "B", "C"), each = 20))
  d$Vol <- rnorm(nrow(d), 50, 10) + rnorm(12, 0, 6)[d$Site]
  d$y <- 1.5 * (d$Loc == "C") + 0.04 * d$Vol +
    rnorm(12, 0, 1.2)[d$Site] + rnorm(nrow(d))
  d
}

test_that("split-plot ANCOVA matches aov for the factor and the covariate", {
  d <- .mixed_cov_example()
  fm <- permanova_mixed_covariate(
    compile_permanova_design(~ Loc + Vol + (1 | Site), d, dist(d$y)),
    permutations = 0
  )
  ref <- summary(aov(y ~ Vol + Loc + Error(Site), d))
  site <- ref[["Error: Site"]][[1L]]
  within <- ref[["Error: Within"]][[1L]]

  expect_equal(fm$table$ss[fm$table$term == "Loc"], site["Loc", "Sum Sq"],
               tolerance = 1e-8)
  expect_equal(fm$table$pseudo_f[fm$table$term == "Loc"], site["Loc", "F value"],
               tolerance = 1e-6)
  expect_equal(fm$table$pseudo_f[fm$table$term == "Vol"], within["Vol", "F value"],
               tolerance = 1e-6)
  # the between-group stratum equals the aov Site-error residual
  stratum <- fm$table[grepl("between-group stratum", fm$table$term), ]
  expect_equal(stratum$df, site["Residuals", "Df"])
  expect_equal(stratum$ms, site["Residuals", "Mean Sq"], tolerance = 1e-8)
})

test_that("covariate-by-factor interaction (slopes) is tested over the residual", {
  d <- .mixed_cov_example()
  fm <- permanova_mixed_covariate(
    compile_permanova_design(~ Loc * Vol + (1 | Site), d, dist(d$y)),
    permutations = 0
  )
  ref <- summary(aov(y ~ Vol * Loc + Error(Site), d))[["Error: Within"]][[1L]]
  expect_equal(fm$table$pseudo_f[fm$table$term == "Loc:Vol"],
               ref["Vol:Loc", "F value"], tolerance = 1e-6)
  # the interaction draws its df from the within residual
  expect_true(fm$table$df[fm$table$term == "Loc:Vol"] == 2L)
})

test_that("permanova() routes a random intercept with a covariate here", {
  d <- .mixed_cov_example()
  fm <- permanova(dist(d$y) ~ Loc + Vol + (1 | Site), d, permutations = 99, seed = 1)
  expect_identical(fm$engine, "balanced split-plot mixed")
  expect_true("Residual" %in% fm$table$Source)
})

test_that("multiple whole-plot factors (no covariate) match the split-plot aov", {
  set.seed(1)
  site <- factor(1:24)
  A <- factor(rep(c("a1", "a2"), each = 12))[site]
  B <- factor(rep(c("b1", "b2", "b3"), 8))[site]
  d <- data.frame(site = rep(site, each = 4))
  d$A <- A[d$site]; d$B <- B[d$site]
  d$y <- c(0, .6)[d$A] + c(0, .4, .9)[d$B] +
    rnorm(24, 0, 1)[d$site] + rnorm(nrow(d))
  fm <- permanova(dist(d$y) ~ A + B + (1 | site), d, permutations = 0)
  expect_identical(fm$engine, "balanced split-plot mixed")
  ref <- summary(aov(y ~ A + B + Error(site), d))[["Error: site"]][[1L]]
  expect_equal(fm$table$`Pseudo-F`[fm$table$Source == "A"], ref["A", "F value"],
               tolerance = 1e-6)
  expect_equal(fm$table$`Pseudo-F`[fm$table$Source == "B"], ref["B", "F value"],
               tolerance = 1e-6)
  # the whole-plot factors are tested over the between-unit stratum
  stratum <- fm$fit$table$ms[grepl("stratum", fm$fit$table$term)]
  expect_equal(stratum, ref["Residuals", "Mean Sq"], tolerance = 1e-8)
})

test_that("a randomized-block factor with no covariate is tested over the residual", {
  set.seed(4)
  nb <- 8
  d <- expand.grid(Treat = factor(c("c", "lo", "hi")), Block = factor(1:nb), rep = 1:3)
  d$y <- c(0, .8, 1.6)[d$Treat] + rnorm(nb, 0, 1)[d$Block] + rnorm(nrow(d))
  fm <- permanova(dist(d$y) ~ Treat + (1 | Block), d, permutations = 0)
  expect_identical(fm$engine, "balanced split-plot mixed")
  ref <- anova(lm(y ~ Block + Treat, d))
  expect_equal(fm$table$`Pseudo-F`[fm$table$Source == "Treat"],
               ref["Treat", "F value"], tolerance = 1e-6)
  expect_equal(fm$fit$table$denominator[fm$fit$table$term == "Treat"], "MS[Residual]")
})

test_that("a between-by-within interaction is routed to the residual", {
  set.seed(5)
  nb <- 8
  d <- expand.grid(T = factor(c("c", "t")), Block = factor(1:nb), rep = 1:3)
  d$A <- factor(c("a1", "a2")[ceiling(as.integer(d$Block) / 4)])  # whole-plot
  d$y <- c(0, .5)[d$A] + c(0, .7)[d$T] + .3 * (d$A == "a2" & d$T == "t") +
    rnorm(nb, 0, .8)[d$Block] + rnorm(nrow(d))
  fm <- permanova(dist(d$y) ~ A * T + (1 | Block), d, permutations = 0)
  ref <- summary(aov(y ~ A * T + Error(Block), d))[["Error: Within"]][[1L]]
  expect_equal(fm$table$`Pseudo-F`[fm$table$Source == "A:T"],
               ref["A:T", "F value"], tolerance = 1e-6)
  # A is whole-plot (unit stratum); T and A:T are within (residual)
  denom <- fm$fit$table$denominator
  expect_equal(denom[fm$fit$table$term == "A:T"], "MS[Residual]")
  expect_equal(denom[fm$fit$table$term == "T"], "MS[Residual]")
})

test_that("a partially confounded factor is refused, not given a wrong denominator", {
  set.seed(2)
  d <- expand.grid(rep = 1:3, Site = factor(1:6))
  # P varies within some sites but not others -> neither whole-plot nor crossed
  d$P <- factor(ifelse(d$Site %in% c("1", "2", "3") & d$rep == 1, "y", "x"))
  d$y <- rnorm(nrow(d))
  expect_error(
    permanova(dist(d$y) ~ P + (1 | Site), d, permutations = 0),
    "partial confounding"
  )
})

test_that("a single whole-plot factor without a covariate stays on permanova_mixed", {
  d <- .mixed_cov_example()
  fm <- permanova(dist(d$y) ~ Loc + (1 | Site), d, permutations = 0)
  expect_identical(fm$engine, "balanced random intercept")
})

test_that("Type I/II/III fixed-effect SS match car::Anova", {
  skip_if_not_installed("car")
  # equal block sizes but an unbalanced (full-rank) A*B cell pattern per block,
  # so the fixed structure is non-orthogonal and the three types differ.
  set.seed(8)
  nb <- 12
  cells <- data.frame(
    A = factor(c("a1", "a1", "a1", "a1", "a2", "a2", "a2"), levels = c("a1", "a2")),
    B = factor(c("b1", "b1", "b2", "b3", "b1", "b2", "b3"), levels = c("b1", "b2", "b3"))
  )
  d <- do.call(rbind, lapply(1:nb, function(b) transform(cells, Block = factor(b))))
  d$Vol <- rnorm(nrow(d), 10, 3)
  d$y <- c(0, 0.7)[d$A] + c(0, 0.5, 1.1)[d$B] + 0.1 * d$Vol +
    rnorm(nb, 0, 0.8)[d$Block] + rnorm(nrow(d))
  getF <- function(tab, term) tab[rownames(tab) == term, "F value"]
  csum <- list(A = contr.sum, B = contr.sum, Block = contr.sum)

  # Type I sequential
  fmI <- permanova(dist(d$y) ~ A * B + Vol + (1 | Block), d, type = "I", permutations = 0)
  refI <- anova(lm(y ~ Block + A + B + A:B + Vol, d))
  expect_equal(fmI$table$`Pseudo-F`[fmI$table$Source == "A"], getF(refI, "A"), tolerance = 1e-6)
  expect_equal(fmI$table$`Pseudo-F`[fmI$table$Source == "A:B"], getF(refI, "A:B"), tolerance = 1e-6)

  # Type II marginal
  fmII <- permanova(dist(d$y) ~ A * B + Vol + (1 | Block), d, type = "II", permutations = 0)
  refII <- car::Anova(lm(y ~ Block + A * B + Vol, d, contrasts = csum), type = 2)
  expect_equal(fmII$table$`Pseudo-F`[fmII$table$Source == "A"], getF(refII, "A"), tolerance = 1e-6)
  expect_equal(fmII$table$`Pseudo-F`[fmII$table$Source == "B"], getF(refII, "B"), tolerance = 1e-6)

  # Type III partial (front door applies contr.sum)
  fmIII <- permanova(dist(d$y) ~ A * B + Vol + (1 | Block), d, type = "III", permutations = 0)
  refIII <- car::Anova(lm(y ~ Block + A * B + Vol, d, contrasts = csum), type = 3)
  expect_equal(fmIII$table$`Pseudo-F`[fmIII$table$Source == "A"], getF(refIII, "A"), tolerance = 1e-6)
  expect_equal(fmIII$table$`Pseudo-F`[fmIII$table$Source == "B"], getF(refIII, "B"), tolerance = 1e-6)
})

test_that("permanova() does not warn that `type` is ignored for this engine", {
  d <- .mixed_cov_example()
  expect_no_warning(
    permanova(dist(d$y) ~ Loc + Vol + (1 | Site), d, type = "II", permutations = 0)
  )
})

test_that("a within-block (sub-plot) factor is tested over the within residual (RCB ANCOVA)", {
  set.seed(4)
  nb <- 8
  d <- expand.grid(Treat = factor(c("ctrl", "lo", "hi")), Block = factor(1:nb), rep = 1:3)
  d$Vol <- rnorm(nrow(d), 20, 5) + rnorm(nb, 0, 3)[d$Block]
  d$y <- c(0, 0.8, 1.6)[d$Treat] + 0.05 * d$Vol +
    rnorm(nb, 0, 1)[d$Block] + rnorm(nrow(d))
  fm <- permanova_mixed_covariate(
    compile_permanova_design(~ Treat + Vol + (1 | Block), d, dist(d$y)),
    permutations = 0
  )
  # matching sequential order: Block, then Treat, then Vol
  ref <- anova(lm(y ~ Block + Treat + Vol, d))
  expect_equal(fm$table$ss[fm$table$term == "Treat"], ref["Treat", "Sum Sq"],
               tolerance = 1e-8)
  expect_equal(fm$table$pseudo_f[fm$table$term == "Treat"], ref["Treat", "F value"],
               tolerance = 1e-6)
  expect_equal(fm$table$pseudo_f[fm$table$term == "Vol"], ref["Vol", "F value"],
               tolerance = 1e-6)
  expect_equal(fm$table$df[fm$table$term == "Residual"], ref["Residuals", "Df"])
  # both within-unit terms are charged to the residual
  expect_true(all(fm$table$denominator[fm$table$term %in% c("Treat", "Vol")] ==
                    "MS[Residual]"))
})

test_that("decomposed covariates are assigned to their actual mixed-model strata", {
  d <- .mixed_cov_example(2)
  d <- decompose_covariate(d, "Vol", by = "Site")
  fm <- permanova(
    dist(d$y) ~ Loc + Vol_between + Vol_within + (1 | Site), d,
    permutations = 0, mc_permutations = 0
  )
  tab <- fm$fit$table
  expect_match(tab$denominator[tab$term == "Vol_between"], "between-group")
  expect_identical(tab$denominator[tab$term == "Vol_within"], "MS[Residual]")
  residual_H <- fm$fit$tests$Vol_within$denominator$projector
  expect_lt(max(abs(residual_H %*% residual_H - residual_H)), 1e-8)
  expect_gte(min(eigen(residual_H, symmetric = TRUE, only.values = TRUE)$values),
             -1e-8)
})

test_that("mixed covariate exact requests enumerate each exchangeability action", {
  set.seed(31)
  d <- expand.grid(rep = 1:2, Site = factor(1:4))
  d$Trt <- factor(rep(c("A", "A", "B", "B"), each = 2))
  d$x <- rep(c(-1, 1), 4)
  d$y <- rnorm(nrow(d))
  fm <- permanova(
    dist(d$y) ~ Trt + x + (1 | Site), d,
    permutations = "all", mc_permutations = 0, seed = 2,
    keep_permutations = TRUE, max_exact = 100000
  )
  expect_equal(fm$fit$tests$Trt$attainable_permutations, 24)
  expect_equal(fm$fit$tests$Trt$permutations_used, 24)
  expect_true(fm$fit$tests$Trt$exact_permutation)
  expect_equal(fm$fit$tests$x$attainable_permutations, 16)
  expect_equal(fm$fit$tests$x$permutations_used, 16)
  expect_true(fm$fit$tests$x$exact_permutation)
})

test_that("a covariate-by-whole-plot interaction names the term that spans it", {
  set.seed(3)
  dat <- do.call(rbind, lapply(1:18, function(i)
    data.frame(Tile = factor(i), Timepoint = c(1, 1, 1, 2, 2))))
  dat$Caging <- factor(rep(c("A", "B", "C"), each = 30))
  dat$y <- rnorm(90) + as.numeric(dat$Caging) + dat$Timepoint
  D <- dist(dat$y)

  # Every tile shares the same mean Timepoint, so the between-unit part of
  # Caging:Timepoint is the Caging indicator scaled -- Caging has nothing left.
  # The refusal is correct; the message has to say which term took the space.
  expect_error(
    permanova(~ Caging * Timepoint + (1 | Tile), dat, distance = D,
              permutations = 0, type = "III"),
    "Caging:Timepoint"
  )
  expect_error(
    permanova(~ Caging * Timepoint + (1 | Tile), dat, distance = D,
              permutations = 0, type = "III"),
    "declare it as a factor"
  )

  # Occasions declared as a factor: the interaction is a within-unit term and
  # the design is an ordinary split-plot.
  dat$Timepoint <- factor(dat$Timepoint)
  fit <- permanova(~ Caging * Timepoint + (1 | Tile), dat, distance = D,
                   permutations = 0, type = "III")
  expect_equal(fit$table$df[1:3], c(2L, 1L, 2L))

  # A covariate that genuinely varies between units stays estimable.
  gradient <- do.call(rbind, lapply(1:18, function(i)
    data.frame(Tile = factor(i), Timepoint = c(1, 1, 1, 2, 2) + (i %% 3))))
  gradient$Caging <- factor(rep(c("A", "B", "C"), each = 30))
  gradient$y <- rnorm(90) + as.numeric(gradient$Caging) + gradient$Timepoint
  expect_s3_class(
    permanova(~ Caging * Timepoint + (1 | Tile), gradient,
              distance = dist(gradient$y), permutations = 0, type = "III"),
    "permanova"
  )
})
