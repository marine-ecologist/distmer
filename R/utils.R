# Shared internal utilities used across engines. Defined once here so the same
# helpers are not redeclared in individual engine files.

#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a

#' @keywords internal
#' @noRd
.stopf <- function(...) stop(sprintf(...), call. = FALSE)

#' @keywords internal
#' @noRd
.warnf <- function(...) warning(sprintf(...), call. = FALSE)
