test_that("mixed formula compiles to a canonical nested design", {
  dat <- expand.grid(
    treatment = factor(c("control", "impact")),
    site = factor(seq_len(3)), block = factor(seq_len(2)), replicate = seq_len(2)
  )
  D <- dist(cbind(seq_len(nrow(dat)), as.numeric(dat$treatment)))
  design <- compile_permanova_design(~ treatment + (1 | site / block), dat, distance = D)

  expect_s3_class(design, "permanova_design")
  expect_equal(design$fixed_effects, "treatment")
  expect_equal(design$random_effects, c("site", "site:block"))
  expect_equal(design$containment_graph$parent, "site")
  expect_equal(design$containment_graph$child, "site:block")
  expect_equal(nrow(design$X), nrow(dat))
  expect_equal(nrow(design$Z[["site:block"]]), nrow(dat))
  expect_true(design$spectral$euclidean)
  expect_equal(design$distance, as.matrix(D))
  expect_true(design$estimability$terms[["treatment"]])
  expect_equal(design$balance$random_balanced, c(site = TRUE, `site:block` = TRUE))
})

test_that("non-Euclidean diagnostics retain negative eigenvalues", {
  D <- matrix(c(0, 1, 1, 1, 0, 3, 1, 3, 0), 3, 3)
  diagnostics <- distance_diagnostics(D)
  expect_false(diagnostics$euclidean)
  expect_gt(diagnostics$negative_count, 0)
  expect_gt(diagnostics$negative_inertia, 0)
  expect_equal(length(diagnostics$eigenvalues), 3)
})

test_that("random slopes require the explicit supported declaration", {
  dat <- data.frame(y = 1:6, treatment = gl(2, 3), site = gl(3, 2))
  expect_error(
    compile_permanova_design(y ~ treatment + (treatment | site), dat),
    "explicit form"
  )

  repeated <- expand.grid(
    replicate = 1:2, condition = factor(c("A", "B")), subject = factor(1:4)
  )
  design <- compile_permanova_design(
    ~ condition + (1 + condition | subject), repeated,
    dist(seq_len(nrow(repeated)))
  )
  expect_equal(design$random_term_table$origin, "slope")
  expect_equal(design$random_term_table$slope, "condition")
  expect_true(design$random_term_table$correlated)
  expect_equal(ncol(design$Z[[1L]]), 8)
})

test_that("explicit random interactions are not mistaken for nesting", {
  dat <- expand.grid(
    replicate = factor(1:2), treatment = factor(c("a", "b")),
    block = factor(1:3), observer = factor(1:2)
  )
  design <- compile_permanova_design(
    ~ treatment + (1 | block) + (1 | observer) + (1 | block:observer),
    dat, dist(seq_len(nrow(dat)))
  )

  expect_equal(design$random_term_table$origin,
               c("main", "main", "interaction"))
  expect_true(is.na(design$random_term_table$parent[[3L]]))
  expect_equal(nrow(design$containment_graph), 0)
  expect_equal(design$interaction_graph$random$term, "block:observer")
  relation <- subset(design$crossing, first == "block" & second == "observer")
  expect_true(relation$crossed)
  expect_equal(relation$relation, "crossed")
})

test_that("slash nesting retains containment provenance", {
  dat <- expand.grid(replicate = factor(1:2), site = factor(1:3), block = factor(1:2))
  design <- compile_permanova_design(
    ~ 1 + (1 | site/block), dat, dist(seq_len(nrow(dat)))
  )
  expect_equal(design$random_term_table$origin, c("main", "nesting"))
  expect_equal(design$random_term_table$parent[[2L]], "site")
  expect_equal(design$containment_graph,
               data.frame(parent = "site", child = "site:block"))
})

test_that("conflicting flattened grouping provenance is rejected", {
  dat <- expand.grid(replicate = factor(1:2), site = factor(1:3), block = factor(1:2))
  expect_error(
    compile_permanova_design(
      ~ 1 + (1 | site/block) + (1 | site:block), dat,
      dist(seq_len(nrow(dat)))
    ),
    "conflicting"
  )
})

test_that("distance row mismatch fails explicitly", {
  dat <- data.frame(group = gl(2, 2))
  expect_error(
    compile_permanova_design(~ group, dat, dist(1:5)),
    "rows must match"
  )
})
