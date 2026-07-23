test_that("capability registry is complete and resolves public implementations", {
  capabilities <- permanova_capabilities()
  required <- c(
    "engine", "companion", "design_envelope", "estimand", "denominator",
    "permutation_contract", "error_control", "explicit_exclusions"
  )
  expect_equal(names(capabilities), required)
  expect_equal(nrow(capabilities), 21)
  expect_equal(anyDuplicated(capabilities$engine), 0L)
  expect_true(all(vapply(
    capabilities$engine, exists, logical(1), mode = "function",
    envir = asNamespace("distmer"), inherits = FALSE
  )))
  companions <- capabilities$companion[nzchar(capabilities$companion)]
  expect_true(all(vapply(
    companions, exists, logical(1), mode = "function",
    envir = asNamespace("distmer"), inherits = FALSE
  )))
  descriptive <- setdiff(required, c("engine", "companion"))
  expect_true(all(vapply(
    capabilities[descriptive], function(value) all(nzchar(value)), logical(1)
  )))
  expect_match(
    capabilities$error_control[capabilities$engine == "permanova_maxT"],
    "weak FWER", fixed = TRUE
  )
  expect_match(
    capabilities$error_control[
      capabilities$engine == "permanova_pairwise_strong"
    ], "strong FWER", fixed = TRUE
  )
  expect_match(
    capabilities$error_control[capabilities$engine == "permanova_mixed_maxT"],
    "weak FWER", fixed = TRUE
  )
})

test_that("human support matrix names every registered implementation", {
  path <- file.path(testthat::test_path(), "..", "..", "inst",
                    "SUPPORT-MATRIX.md")
  if (!file.exists(path)) {
    path <- system.file("SUPPORT-MATRIX.md", package = "distmer")
  }
  expect_true(nzchar(path) && file.exists(path))
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_true(all(vapply(
    permanova_capabilities()$engine,
    function(engine) grepl(paste0("`", engine, "`"), text, fixed = TRUE),
    logical(1)
  )))
  expect_match(text, "runtime validators remain authoritative", fixed = TRUE)
  expect_match(text, "Strong control is confined", fixed = TRUE)
})
