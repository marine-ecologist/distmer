root <- normalizePath(getwd(), mustWork = TRUE)
description <- read.dcf(file.path(root, "DESCRIPTION"))
version <- unname(description[1L, "Version"])
news <- readLines(file.path(root, "NEWS.md"), warn = FALSE)
news_version <- sub("^# distmer ", "", news[1L])

namespace_lines <- readLines(file.path(root, "NAMESPACE"), warn = FALSE)
export_lines <- namespace_lines[startsWith(namespace_lines, "export(")]
exports <- substring(export_lines, 8L, nchar(export_lines) - 1L)
s3_lines <- namespace_lines[startsWith(namespace_lines, "S3method(")]

rd_files <- list.files(
  file.path(root, "man"), pattern = "[.]Rd$", full.names = TRUE
)
rd <- unlist(lapply(rd_files, readLines, warn = FALSE), use.names = FALSE)
alias_lines <- grep("^\\\\alias\\{", rd, value = TRUE)
aliases <- sub("^\\\\alias\\{([^}]*)\\}$", "\\1", alias_lines)

vignettes <- list.files(
  file.path(root, "vignettes"), pattern = "[.]Rmd$", full.names = TRUE
)
ignore <- readLines(file.path(root, ".Rbuildignore"), warn = FALSE)
authors <- unname(description[1L, "Authors@R"])
license <- readLines(file.path(root, "LICENSE"), warn = FALSE)
support_matrix <- file.path(root, "inst", "SUPPORT-MATRIX.md")
release_handoff <- file.path(root, "inst", "RELEASE-HANDOFF.md")
version_policy <- file.path(root, "inst", "VERSIONING.md")

structural <- c(
  `DESCRIPTION/NEWS version parity` = identical(version, news_version),
  `all exports have Rd aliases` = length(setdiff(exports, aliases)) == 0L,
  `public exports are unique` = !anyDuplicated(exports),
  `S3 registrations are unique` = !anyDuplicated(s3_lines),
  `vignette series contains 7 documents` = length(vignettes) == 7L,
  `validated support matrix is present` = file.exists(support_matrix),
  `development release handoff is present` = file.exists(release_handoff),
  `semantic version policy is present` = file.exists(version_policy),
  `build ignores check and tar artifacts` =
    any(grepl("Rcheck", ignore, fixed = TRUE)) &&
    any(grepl("tar", ignore, fixed = TRUE)),
  `build ignores Finder metadata` =
    any(grepl("DS_Store", ignore, fixed = TRUE))
)

metadata_blockers <- c(
  if (grepl("maintainer@example.org", authors, fixed = TRUE)) {
    "Authors@R still contains maintainer@example.org"
  },
  if (any(grepl("Project Maintainer", license, fixed = TRUE))) {
    "LICENSE still contains the Project Maintainer placeholder"
  },
  if (!"URL" %in% colnames(description)) "DESCRIPTION has no project URL",
  if (!"BugReports" %in% colnames(description)) {
    "DESCRIPTION has no BugReports URL"
  }
)

cat("distmer release audit\n")
cat("Version:", version, "\n")
for (name in names(structural)) {
  cat(if (structural[[name]]) "PASS" else "FAIL", "-", name, "\n")
}
cat("Public exports:", length(exports), "\n")
cat("S3 registrations:", length(s3_lines), "\n")
cat("Vignettes:", length(vignettes), "\n")
if (length(metadata_blockers)) {
  cat("PUBLICATION BLOCKED:\n")
  cat(paste0("- ", metadata_blockers, collapse = "\n"), "\n")
} else {
  cat("Publication metadata: PASS\n")
}

if (!all(structural)) quit(status = 1L)
if (identical(tolower(Sys.getenv("STRICT_RELEASE_METADATA")), "true") &&
    length(metadata_blockers)) quit(status = 2L)
