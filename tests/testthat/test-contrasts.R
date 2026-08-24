# Contrasts on a fitted permanova() model (the fit-based API). The correctness
# anchor is the deterministic identity: method = "pooled" pseudo-F equals the
# classical constrained-scalar contrast F. method = "isolated" runs the
# validated two-group engine on the involved levels.

.contrast_example <- function() {
  data.frame(
    group = factor(c("A", "A", "B", "B", "C"), levels = c("A", "B", "C")),
    y = c(-1.2, 0.2, 0.4, 1.8, 2.5)
  )
}

.direct_contrast_f <- function(y, group, coefficients) {
  means <- tapply(y, group, mean)
  counts <- table(group)
  contrast_ss <- sum(coefficients * means) ^ 2 / sum(coefficients ^ 2 / counts)
  full <- lm(y ~ group)
  contrast_ss / (sum(residuals(full) ^ 2) / df.residual(full))
}

test_that("pooled contrast F equals the scalar contrast formula", {
  dat <- .contrast_example()
  fit <- permanova(dist(dat$y) ~ group, dat, permutations = "all")
  L <- rbind(`A - B` = c(A = 1, B = -1, C = 0),
             `AB - C` = c(A = 0.5, B = 0.5, C = -1))
  ct <- suppressWarnings(
    permanova_contrasts(fit, type = L, method = "pooled", adjust = "none")
  )
  direct <- apply(L, 1L, function(cf) .direct_contrast_f(dat$y, dat$group, cf))
  expect_equal(ct$table$pseudo_F, unname(direct), tolerance = 1e-8)
  expect_s3_class(ct, "permanova_contrasts")
})

test_that("isolated contrasts run the two-group engine on the involved levels", {
  dat <- .contrast_example()
  fit <- permanova(dist(dat$y) ~ group, dat, permutations = "all")
  ct <- suppressWarnings(permanova_contrasts(
    fit, type = rbind(`A - B` = c(A = 1, B = -1, C = 0)), method = "isolated"
  ))
  # a two-level contrast is a valid two-group permutation test
  expect_true(is.finite(ct$table$p_perm[[1L]]))
  # it equals the two-group permanova restricted to A and B
  sub <- dat[dat$group %in% c("A", "B"), ]
  sub$group <- droplevels(sub$group)
  two <- suppressWarnings(permanova(dist(sub$y) ~ group, sub, permutations = "all"))
  expect_equal(ct$table$pseudo_F[[1L]], two$table$`Pseudo-F`[[1L]], tolerance = 1e-8)
})

test_that("the pairwise family is built and multiplicity-adjusted", {
  dat <- .contrast_example()
  fit <- permanova(dist(dat$y) ~ group, dat, permutations = "all")
  ct <- suppressWarnings(
    permanova_contrasts(fit, type = "pairwise", adjust = "bonferroni")
  )
  expect_equal(ct$table$contrast, c("A - B", "A - C", "B - C"))
  estimable <- !is.na(ct$table$p_perm)
  expect_equal(ct$table$p_adj[estimable],
               p.adjust(ct$table$p_perm[estimable], method = "bonferroni"))
})

test_that("contrast_matrix builds named families that sum to zero", {
  levs <- c("A", "B", "C")
  pw <- contrast_matrix("pairwise", levs)
  expect_equal(rownames(pw), c("A - B", "A - C", "B - C"))
  expect_true(all(abs(rowSums(pw)) < 1e-8))
  tvc <- contrast_matrix("trt.vs.ctrl", levs, ref = "A")
  expect_true(all(abs(rowSums(tvc)) < 1e-8))
})

test_that("permanova_contrasts requires a fitted permanova object", {
  dat <- .contrast_example()
  design <- compile_permanova_design(~ group, dat, dist(dat$y))
  expect_error(permanova_contrasts(design, type = "pairwise"),
               "permanova\\(\\) object")
})

test_that("non-Euclidean dissimilarities keep empirical isolated inference", {
  dat <- .contrast_example()
  D <- as.matrix(dist(cbind(dat$y, seq_len(nrow(dat)))))
  D[1L, 2L] <- D[2L, 1L] <- 5 * max(D)
  fit <- permanova(as.dist(D) ~ group, dat, permutations = "all")
  expect_false(fit$design$spectral$euclidean)
  ct <- suppressWarnings(permanova_contrasts(fit, type = "pairwise"))
  expect_true(any(is.finite(ct$table$p_perm)))
})

test_that("isolated contrasts preserve the fitted random-intercept design", {
  set.seed(16)
  d<-expand.grid(rep=1:3,site_id=1:3,trt=factor(c("A","B","C")))
  d$site<-interaction(d$trt,d$site_id,drop=TRUE)
  d$y<-rnorm(nlevels(d$site),sd=2)[d$site]+rnorm(nrow(d),sd=.2)
  fit<-permanova(dist(d$y)~trt+(1|site),d,permutations=99,seed=2)
  ct<-suppressWarnings(permanova_contrasts(fit,type="pairwise",permutations=99,seed=2))
  s<-droplevels(d[d$trt%in%c("A","B"),])
  ref<-permanova(dist(s$y)~trt+(1|site),s,permutations=99,seed=2)
  row<-ct$table[ct$table$contrast=="A - B",]
  expect_equal(row$pseudo_F,ref$table$`Pseudo-F`[1],tolerance=1e-10)
  expect_equal(row$p_perm,ref$table$`P(perm)`[1])
  expect_equal(row$df_den,.permanova_tests(ref$fit)$trt$denominator$df)
})


# ---------------------------------------------------------------------------
# within-cell dispatch: permanova_contrasts(within = ) and the multi-name
# `factor` shorthand both route to permanova_contrasts_within().
# ---------------------------------------------------------------------------

.within_example <- function(seed = 4) {
  set.seed(seed)
  dat <- expand.grid(
    replicate = 1:3,
    block = factor(1:4),
    depth = factor(c("lo", "hi")),
    region = factor(c("north", "south"))
  )
  dat$block <- factor(paste(dat$region, dat$block))
  dat$y <- rnorm(nrow(dat)) + as.numeric(dat$depth) +
    as.numeric(dat$block) / 4
  dat
}

test_that("a multi-name `factor` is the shorthand for `within`", {
  dat <- .within_example()
  fit <- permanova(dist(dat$y) ~ depth * region + (1 | block), dat,
                   permutations = 0)

  expect_message(
    shorthand <- permanova_contrasts(fit, type = "pairwise",
                                     factor = c("depth", "region")),
    "Contrasting 'depth' within region"
  )
  explicit <- permanova_contrasts(fit, type = "pairwise", factor = "depth",
                                  within = "region")
  expect_equal(as.data.frame(shorthand), as.data.frame(explicit))
  expect_s3_class(shorthand, "permanova_contrasts_within")

  # one row per pair per cell, and the cell columns are carried through
  expect_true(all(c("region", "contrast", "t", "denominator_df") %in%
                    names(shorthand)))
  expect_equal(nrow(shorthand), nlevels(dat$region))
  expect_setequal(as.character(shorthand$region), levels(dat$region))
})

test_that("the within path refuses families and denominators it cannot form", {
  dat <- .within_example()
  fit <- permanova(dist(dat$y) ~ depth * region + (1 | block), dat,
                   permutations = 0)

  expect_error(
    permanova_contrasts(fit, type = "poly", factor = "depth",
                        within = "region"),
    'only defined for type = "pairwise"'
  )
  expect_error(
    permanova_contrasts(fit, type = "pairwise", method = "pooled",
                        factor = "depth", within = "region"),
    "no pooled form"
  )
  # an unknown name is reported, not silently announced as a within-cell run
  expect_error(
    suppressMessages(
      permanova_contrasts(fit, type = "pairwise", factor = c("depth", "nope"))
    ),
    "nope"
  )
})

test_that("printing a within-cell table survives losing its cell columns", {
  dat <- .within_example()
  fit <- permanova(dist(dat$y) ~ depth * region + (1 | block), dat,
                   permutations = 0)
  ct <- permanova_contrasts(fit, type = "pairwise", factor = "depth",
                            within = "region")
  expect_output(print(ct), "Pairwise contrasts of 'depth' within region")
  # subsetting away the grouping columns must fall back, not fail
  expect_output(print(ct[, c("contrast", "t")]), "contrast")
})
