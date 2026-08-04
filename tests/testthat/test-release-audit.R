.package_root <- function() {
  candidate <- normalizePath(
    file.path(testthat::test_path(), "..", ".."), mustWork = TRUE
  )
  if (file.exists(file.path(candidate, "DESCRIPTION"))) candidate else
    system.file(package = "distmer")
}

.package_file <- function(...) {
  path <- file.path(.package_root(), ...)
  if (file.exists(path)) path else system.file(..., package = "distmer")
}

test_that("release version and NEWS heading agree", {
  root <- .package_root()
  description <- read.dcf(.package_file("DESCRIPTION"))
  news <- readLines(.package_file("NEWS.md"), warn = FALSE)
  expect_equal(unname(description[1L, "Version"]),
               sub("^# distmer ", "", news[1L]))
})

test_that("every namespace export has a documentation alias", {
  root <- .package_root()
  namespace <- readLines(.package_file("NAMESPACE"), warn = FALSE)
  export_lines <- namespace[startsWith(namespace, "export(")]
  exports <- substring(export_lines, 8L, nchar(export_lines) - 1L)
  rd_paths <- list.files(file.path(root, "man"), pattern = "[.]Rd$",
                         full.names = TRUE)
  if (length(rd_paths)) {
    rd <- unlist(lapply(rd_paths, readLines, warn = FALSE), use.names = FALSE)
    alias_lines <- grep("^\\\\alias\\{", rd, value = TRUE)
    aliases <- sub("^\\\\alias\\{([^}]*)\\}$", "\\1", alias_lines)
  } else {
    index <- readLines(.package_file("help", "AnIndex"), warn = FALSE)
    aliases <- sub("\\t.*$", "", index)
  }

  expect_setequal(setdiff(exports, aliases), character())
  expect_equal(anyDuplicated(exports), 0L)
})

test_that("registered S3 methods resolve", {
  root <- .package_root()
  namespace <- readLines(.package_file("NAMESPACE"), warn = FALSE)
  lines <- namespace[startsWith(namespace, "S3method(")]
  contents <- substring(lines, 10L, nchar(lines) - 1L)
  entries <- strsplit(contents, ",", fixed = TRUE)
  # A generic may be namespace-qualified (e.g. `base::print` for a method on a
  # base generic); the method object itself is `generic.class`, so strip the
  # `pkg::` qualifier before resolving.
  methods <- vapply(entries, function(entry) {
    entry <- trimws(entry)
    entry[[1L]] <- sub("^.*::", "", entry[[1L]])
    paste(entry, collapse = ".")
  }, character(1))
  namespace_environment <- asNamespace("distmer")
  expect_true(all(vapply(
    methods, exists, logical(1), envir = namespace_environment,
    inherits = FALSE
  )))
  expect_equal(anyDuplicated(lines), 0L)
})

test_that("release artifacts and platform metadata are excluded", {
  root <- .package_root()
  ignore_path <- file.path(root, ".Rbuildignore")
  if (file.exists(ignore_path)) {
    ignore <- readLines(ignore_path, warn = FALSE)
    expect_true(any(grepl("Rcheck", ignore, fixed = TRUE)))
    expect_true(any(grepl("tar", ignore, fixed = TRUE)))
    expect_true(any(grepl("DS_Store", ignore, fixed = TRUE)))
  } else {
    succeed("Source-only build-ignore rules were applied before installation.")
  }
  source_vignettes <- file.path(root, "vignettes")
  doc_path <- if (dir.exists(source_vignettes)) source_vignettes else
    .package_file("doc")
  expect_length(list.files(doc_path, "[.]Rmd$"), 7)
})

test_that("publication metadata is present", {
  root <- .package_root()
  readiness_path <- file.path(root, "inst", "RELEASE-READINESS.md")
  if (!file.exists(readiness_path)) {
    readiness_path <- .package_file("RELEASE-READINESS.md")
  }
  readiness <- paste(readLines(readiness_path, warn = FALSE), collapse = "\n")
  expect_match(readiness, "BugReports", fixed = TRUE)
  expect_match(readiness, "project `URL`", fixed = TRUE)
  expect_match(readiness, "strict metadata audit", fixed = TRUE)
})

test_that("development handoff preserves the owner publication boundary", {
  root <- .package_root()
  handoff_path <- file.path(root, "inst", "RELEASE-HANDOFF.md")
  if (!file.exists(handoff_path)) {
    handoff_path <- .package_file("RELEASE-HANDOFF.md")
  }
  handoff <- paste(readLines(handoff_path, warn = FALSE), collapse = "\n")
  description <- read.dcf(.package_file("DESCRIPTION"))
  version <- unname(description[1L, "Version"])
  expect_match(handoff, paste(version, "development tree"), fixed = TRUE)
  expect_match(handoff, "Only the repository owner", fixed = TRUE)
  expect_match(handoff, "STRICT_RELEASE_METADATA=true", fixed = TRUE)
})

test_that("bounded benchmark smoke completes", {
  timing <- benchmark_permanova(
    group_counts = 4L, replicates = 2L, permutations = 9L, seed = 20260719
  )
  expect_equal(timing$observations, 8)
  expect_true(is.finite(timing$elapsed_seconds))
  expect_lt(timing$elapsed_seconds, 30)
})
