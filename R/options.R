#' Options that change how distmer reports
#'
#' `distmer` never changes an inferential result through an option: every one
#' below controls reporting only. Set them with [options()], for example
#' `options(distmer.resolution_alpha = 0.01)`.
#'
#' @section Permutation resolution:
#'
#' `distmer.resolution_alpha` (default `0.05`) is the level a permutation
#' distribution's *resolution* is judged against. A restricted permutation
#' scheme can attain many arrangements while producing only a handful of
#' distinct statistics -- permuting labels inside a treatment group leaves the
#' statistic untouched -- so the smallest p-value a design can return may sit
#' well above zero. Two groups of two attain 24 arrangements but only three
#' distinct pseudo-F values, and can never return a p-value below a third.
#'
#' Every test specification records `distinct_statistics` and
#' `minimum_attainable_p` regardless of this option. The option decides only
#' when a design that cannot reach the chosen level says so in its warnings:
#' lower it to `0.01` to be told only about designs that cannot reach one
#' percent, or raise it to `0.1` to be warned more readily.
#'
#' @section Progress reporting:
#'
#' `distmer.progress` (default: whether the session is interactive) turns the
#' permutation progress bar on or off. `distmer.progress_min` (default `200`)
#' is the number of permutations below which a run is too short to be worth
#' reporting. Both are ignored outside an interactive session unless
#' `distmer.progress` is set explicitly, so scripts, tests, vignettes, and
#' `R CMD check` stay quiet.
#'
#' @name distmer-options
#' @seealso [permanova()]
NULL
