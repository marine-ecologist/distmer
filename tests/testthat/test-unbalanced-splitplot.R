# Unbalanced multi-factor split-plot: whole-plot factors are tested on equally
# weighted unit centroids (matching an lm on island means), and within-unit
# factors over the within-unit residual (matching a unit-blocked lm). For a
# balanced design it reduces to permanova_mixed_covariate.

.splitplot_example <- function(seed = 3) {
  set.seed(seed)
  cells <- list(c("A", "unpop"), c("A", "pop"), c("B", "unpop"), c("B", "pop"))
  islands_per_cell <- c(2, 3, 2, 3)
  rows <- list(); id <- 0L
  for (c in seq_along(cells)) for (i in seq_len(islands_per_cell[[c]])) {
    id <- id + 1L
    for (de in c("shallow", "deep")) {
      nrep <- sample(2:4, 1L)
      rows[[length(rows) + 1L]] <- data.frame(
        EC = cells[[c]][[1L]], PS = cells[[c]][[2L]],
        ISLAND = paste0("i", id), DE = de, r = seq_len(nrep)
      )
    }
  }
  d <- do.call(rbind, rows)
  for (v in c("EC", "PS", "DE", "ISLAND")) d[[v]] <- factor(d[[v]])
  d$y <- 0.8 * (d$EC == "B") + 0.5 * (d$DE == "deep") +
    0.4 * (d$DE == "deep" & d$EC == "B") +
    rnorm(nlevels(d$ISLAND), 0, 1.2)[d$ISLAND] + rnorm(nrow(d))
  d
}

test_that("whole-plot factors match an lm on equal-weight unit centroids", {
  d <- .splitplot_example()
  fm <- permanova_unbalanced_splitplot(
    compile_permanova_design(~ EC * PS + DE + (1 | ISLAND), d, dist(d$y)),
    permutations = 0
  )
  island_means <- aggregate(y ~ ISLAND + EC + PS, d, mean)
  ref <- anova(lm(y ~ EC * PS, island_means))
  for (term in c("EC", "PS", "EC:PS")) {
    expect_equal(fm$table$pseudo_f[fm$table$term == term],
                 ref[term, "F value"], tolerance = 1e-6)
  }
  # tested over the between-unit (island) stratum
  expect_true(all(fm$table$denominator[fm$table$term %in% c("EC", "PS", "EC:PS")] ==
                    "MS[ISLAND stratum]"))
})

test_that("within-unit factors and interactions match a unit-blocked lm", {
  d <- .splitplot_example()
  fm <- permanova_unbalanced_splitplot(
    compile_permanova_design(~ EC * PS * DE + (1 | ISLAND), d, dist(d$y)),
    permutations = 0
  )
  ref <- anova(lm(y ~ ISLAND + DE + DE:EC + DE:PS + DE:EC:PS, d))
  pairs <- c(DE = "DE", `EC:DE` = "DE:EC", `PS:DE` = "DE:PS", `EC:PS:DE` = "DE:EC:PS")
  for (nm in names(pairs)) {
    expect_equal(fm$table$pseudo_f[fm$table$term == nm],
                 ref[pairs[[nm]], "F value"], tolerance = 1e-6)
  }
  expect_true(all(fm$table$denominator[fm$table$term %in% names(pairs)] ==
                    "MS[Residual]"))
})

test_that("a balanced design reduces to permanova_mixed_covariate", {
  set.seed(1)
  b <- expand.grid(r = 1:3, DE = factor(c("s", "d")), isl = 1:2,
                   PS = factor(c("u", "p")), EC = factor(c("A", "B")))
  b$ISLAND <- interaction(b$EC, b$PS, b$isl, drop = TRUE)
  b$y <- 0.8 * (b$EC == "B") + 0.5 * (b$DE == "d") +
    rnorm(nlevels(b$ISLAND), 0, 1.2)[b$ISLAND] + rnorm(nrow(b))
  cm <- compile_permanova_design(~ EC * PS * DE + (1 | ISLAND), b, dist(b$y))
  unbal <- permanova_unbalanced_splitplot(cm, permutations = 0)$table
  bal <- permanova_mixed_covariate(cm, permutations = 0)$table
  for (term in bal$term[!is.na(bal$pseudo_f)]) {
    expect_equal(unbal$pseudo_f[unbal$term == term],
                 bal$pseudo_f[bal$term == term], tolerance = 1e-8)
  }
})

test_that("a single whole-plot factor matches permanova_unequal_random_intercept", {
  d <- .splitplot_example()
  cm <- compile_permanova_design(~ EC + (1 | ISLAND), d, dist(d$y))
  expect_equal(
    permanova_unbalanced_splitplot(cm, permutations = 0)$table$pseudo_f[1],
    permanova_unequal_random_intercept(cm, permutations = 0)$table$pseudo_f[1],
    tolerance = 1e-8
  )
})

test_that("permanova() routes an unbalanced multi-factor mixed design here", {
  d <- .splitplot_example()
  fm <- permanova(dist(d$y) ~ EC * PS * DE + (1 | ISLAND), d,
                  permutations = 99, seed = 1)
  expect_identical(fm$engine, "unbalanced split-plot mixed")
  a <- permanova(dist(d$y) ~ EC * PS + DE + (1 | ISLAND), d, permutations = 99, seed = 7)
  b <- permanova(dist(d$y) ~ EC * PS + DE + (1 | ISLAND), d, permutations = 99, seed = 7)
  expect_equal(a$table$`P(perm)`, b$table$`P(perm)`)
})

test_that("a sub-plot factor missing from some units still runs when estimable", {
  d <- .splitplot_example()
  # island i1 loses its 'deep' samples -> not fully crossed, but DE is still
  # estimable within the other islands.
  d2 <- droplevels(d[!(d$ISLAND == "i1" & d$DE == "deep"), ])
  fm <- permanova_unbalanced_splitplot(
    compile_permanova_design(~ EC * PS + DE + (1 | ISLAND), d2, dist(d2$y)),
    permutations = 0
  )
  expect_equal(fm$table$denominator[fm$table$term == "DE"], "MS[Residual]")
  ref <- anova(lm(y ~ ISLAND + DE, d2))
  expect_equal(fm$table$pseudo_f[fm$table$term == "DE"],
               ref["DE", "F value"], tolerance = 1e-6)
})

test_that("a genuinely confounded factor (within rank below its rank) is refused", {
  set.seed(9)
  isl <- rep(1:6, each = 4)
  Z <- character(length(isl))
  for (i in 1:4) Z[isl == i] <- rep(c("z1", "z2"), 2)  # islands 1-4 vary z1/z2
  Z[isl %in% c(5, 6)] <- "z3"                            # islands 5-6 are all z3
  d <- data.frame(ISLAND = factor(isl),
                  EC = factor(ifelse(isl <= 3, "A", "B")), Z = factor(Z))
  d$y <- rnorm(length(isl)) + rnorm(6, 0, 1)[d$ISLAND]
  expect_error(
    permanova_unbalanced_splitplot(
      compile_permanova_design(~ EC + Z + (1 | ISLAND), d, dist(d$y)), permutations = 0
    ), "partial confounding"
  )
})

test_that("exact split-plot actions and size strata are enforced", {
  set.seed(13)
  d <- expand.grid(DE=factor(c("s","d")), ISLAND=factor(1:4))
  d$TRT <- factor(rep(c("A","A","B","B"), each=2)); d$y<-rnorm(nrow(d))
  fm <- permanova_unbalanced_splitplot(compile_permanova_design(
    ~ TRT+DE+(1|ISLAND), d, dist(d$y)), permutations="all", max_exact=100,
    keep_permutations=TRUE)
  expect_length(fm$tests$TRT$permutation_statistics,24)
  expect_length(fm$tests$DE$permutation_statistics,16)
  sizes <- rep(2:4,each=2)
  rows <- lapply(seq_along(sizes), function(i) data.frame(
    ISLAND=factor(i), TRT=factor(ifelse(i%%2,"A","B")), y=rnorm(sizes[i])))
  z<-droplevels(do.call(rbind,rows))
  sz<-permanova_unbalanced_splitplot(compile_permanova_design(
    ~TRT+(1|ISLAND),z,dist(z$y)),permutations="all",max_exact=100,keep_permutations=TRUE)
  expect_length(sz$tests$TRT$permutation_statistics,8)
})
