#' Render a table with kableExtra styling
#'
#' Convenience wrapper around `kableExtra::kbl()` and `kableExtra::kable_styling()`
#' for quick, consistently formatted HTML tables in R Markdown/Quarto output.
#'
#' @param x A data frame, matrix, or other object accepted by `kbl()`.
#' @param digits Number of decimal places for numeric columns. Default 3.
#' @param ... Additional arguments passed to `kbl()`.
#' @return A styled `kableExtra` table object.
#' @importFrom kableExtra kbl kable_styling
#' @examples
#' tab(head(mtcars))
#' @export
tab <- function(x, digits = 3, ...) {
  x |>
    kbl(digits = digits, ...) |>
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"))
}
