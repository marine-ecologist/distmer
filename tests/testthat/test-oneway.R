test_that("direct Anderson equations match projection formulation", {
  y <- matrix(
    c(
      0, 0,
      0, 1,
      1, 0,
      4, 4,
      4, 5,
      5, 4
    ),
    ncol = 2,
    byrow = TRUE
  )
  group <- factor(rep(c("A", "B"), each = 3))
  D <- as.matrix(dist(y))

  direct <- permanova_oneway(D, group, permutations = 0, engine = "pairwise")
  projected <- permanova_oneway(D, group, permutations = 0, engine = "projection")

  expect_equal(direct$table$ss, projected$table$ss, tolerance = 1e-10)
  expect_equal(direct$table$ms, projected$table$ms, tolerance = 1e-10)
  expect_equal(direct$table$pseudo_f, projected$table$pseudo_f, tolerance = 1e-10)
})

test_that("Euclidean one-way pseudo-F equals ordinary ANOVA F", {
  y <- c(1.1, 1.5, 1.7, 4.2, 4.8, 5.1)
  group <- factor(rep(c("A", "B"), each = 3))
  D <- as.matrix(dist(y))

  result <- permanova_oneway(D, group, permutations = 0, engine = "projection")
  reference <- summary(aov(y ~ group))[[1]][["F value"]][1]

  expect_equal(result$table$pseudo_f[1], unname(reference), tolerance = 1e-10)
})

test_that("pairwise engine rejects unbalanced replication", {
  y <- matrix(seq_len(10), ncol = 2)
  group <- factor(c("A", "A", "B", "B", "B"))

  expect_error(
    permanova_oneway(dist(y), group, permutations = 0, engine = "pairwise"),
    "balanced"
  )
})
