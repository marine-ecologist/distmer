# Asymmetric nested mixed design (PRIMER ST5 structure): whole-plot fixed over
# IS(unit), sub-plot fixed over subplot x IS(unit), both randoms over residual.
# The univariate anchor is aov(y ~ DE*EC*PS + Error(ISLAND/DE)).

.asym_example <- function(seed = 4) {
  set.seed(seed)
  d <- expand.grid(rep = 1:2, DE = factor(c("s", "m", "d")), isl = 1:2,
                   PS = factor(c("u", "p")), EC = factor(c("A", "B", "C")))
  d$ISLAND <- interaction(d$EC, d$PS, d$isl, drop = TRUE)
  ie <- rnorm(nlevels(d$ISLAND), 0, 1.0)[d$ISLAND]
  cell <- interaction(d$ISLAND, d$DE, drop = TRUE)
  d$y <- 0.7 * (d$EC == "C") + 0.5 * (d$DE == "d") +
    0.4 * (d$DE == "d" & d$PS == "p") + ie +
    rnorm(nlevels(cell), 0, 0.7)[cell] + rnorm(nrow(d), 0, 0.6)
  d
}

test_that("all fixed terms and strata match aov with an ISLAND/DE error structure", {
  d <- .asym_example()
  fm <- permanova_asymmetric_mixed(
    compile_permanova_design(~ DE * EC * PS + (1 | ISLAND) + (1 | ISLAND:DE),
                             d, dist(d$y)), permutations = 0
  )
  a <- summary(aov(y ~ DE * EC * PS + Error(ISLAND / DE), d))
  isl <- a[["Error: ISLAND"]][[1L]]
  dei <- a[["Error: ISLAND:DE"]][[1L]]
  rownames(isl) <- trimws(rownames(isl))
  rownames(dei) <- trimws(rownames(dei))
  fval <- function(term) fm$table$pseudo_f[fm$table$term == term]
  # whole-plot terms over IS(unit)
  for (term in c("EC", "PS", "EC:PS")) {
    expect_equal(fval(term), isl[term, "F value"], tolerance = 1e-6)
    expect_equal(fm$table$denominator[fm$table$term == term], "MS[IS(ISLAND)]")
  }
  # sub-plot terms over subplot x IS(unit)
  for (term in c("DE", "DE:EC", "DE:PS", "DE:EC:PS")) {
    expect_equal(fval(term), dei[term, "F value"], tolerance = 1e-6)
    expect_equal(fm$table$denominator[fm$table$term == term], "MS[DE x IS(ISLAND)]")
  }
  # the two random strata sit over the residual with the aov error df
  expect_equal(fm$table$df[fm$table$term == "IS(ISLAND)"], isl["Residuals", "Df"])
  expect_equal(fm$table$df[fm$table$term == "DE x IS(ISLAND)"], dei["Residuals", "Df"])
  expect_true(all(fm$table$denominator[grepl("IS\\(ISLAND\\)", fm$table$term)] ==
                    "MS[Residual]"))
})

test_that("permanova() routes a unit + unit:subplot design here, reproducibly", {
  d <- .asym_example()
  fm <- permanova(dist(d$y) ~ DE * EC * PS + (1 | ISLAND) + (1 | ISLAND:DE),
                  d, permutations = 199, seed = 1)
  expect_identical(fm$engine, "asymmetric nested mixed")
  b <- permanova(dist(d$y) ~ DE * EC * PS + (1 | ISLAND) + (1 | ISLAND:DE),
                 d, permutations = 199, seed = 1)
  expect_equal(fm$table$`P(perm)`, b$table$`P(perm)`)
})

test_that("degrees of freedom partition every observation", {
  d <- .asym_example()
  fm <- permanova_asymmetric_mixed(
    compile_permanova_design(~ DE * EC * PS + (1 | ISLAND) + (1 | ISLAND:DE),
                             d, dist(d$y)), permutations = 0
  )
  expect_equal(sum(fm$table$df), nrow(d) - 1L)
})

test_that("exact asymmetric requests enumerate and incomplete topologies run", {
  set.seed(11)
  d <- expand.grid(rep=1:2, DE=factor(c("s","d")), ISLAND=factor(1:4))
  d$TRT <- factor(rep(c("A","A","B","B"), each=4)); d$y <- rnorm(nrow(d))
  whole <- permanova_asymmetric_mixed(compile_permanova_design(
    ~ TRT + (1|ISLAND) + (1|ISLAND:DE), d, dist(d$y)),
    permutations="all", max_exact=1000, keep_permutations=TRUE)
  expect_true(whole$tests$TRT$exact_permutation)
  expect_length(whole$tests$TRT$permutation_statistics, 24)
  w <- expand.grid(rep=1:2, DE=factor(c("s","d")), ISLAND=factor(1:2)); w$y<-rnorm(nrow(w))
  within <- permanova_asymmetric_mixed(compile_permanova_design(
    ~ DE + (1|ISLAND) + (1|ISLAND:DE), w, dist(w$y)),
    permutations="all", max_exact=1000, keep_permutations=TRUE)
  expect_length(within$tests$DE$permutation_statistics, 576)
  u <- .asym_example(); first <- levels(u$ISLAND)[1]
  u <- droplevels(u[!(u$ISLAND==first & u$DE==levels(u$DE)[1]),])
  expect_no_error(permanova_asymmetric_mixed(compile_permanova_design(
    ~ DE*EC*PS + (1|ISLAND) + (1|ISLAND:DE), u, dist(u$y)), permutations=19, seed=4))
})
