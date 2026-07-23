test_that("sequential Euclidean sums of squares equal aov for a balanced fixed model", {
  dat <- expand.grid(
    A = factor(c("a1", "a2")),
    B = factor(c("b1", "b2", "b3")),
    replicate = seq_len(3)
  )
  dat$y <- with(dat,
    2 * as.numeric(A) +
      0.5 * as.numeric(B) +
      as.numeric(A) * as.numeric(B) * 0.2 +
      replicate * 0.03
  )

  result <- permanova_fixed(
    distance = dist(dat$y),
    formula = ~ A * B,
    data = dat,
    permutations = 0
  )

  reference <- summary(aov(y ~ A * B, data = dat))[[1]]
  result_terms <- result$table[result$table$term %in% c("A", "B", "A:B", "Residual"), ]

  expect_equal(result_terms$ss, reference[["Sum Sq"]], tolerance = 1e-8)
  expect_equal(result_terms$df, reference[["Df"]], tolerance = 1e-8)
  expect_equal(result_terms$pseudo_f[1:3], reference[["F value"]][1:3], tolerance = 1e-8)
})
