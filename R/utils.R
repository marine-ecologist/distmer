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

# Permutation p-value plus the resolution of the distribution it came from.
#
# A restricted permutation scheme can attain many arrangements yet produce only
# a handful of distinct statistics: permuting labels inside a treatment group
# leaves the statistic untouched, so the achievable p-values are coarse. Two
# groups of two attain 24 arrangements but only three distinct pseudo-F values,
# and the smallest p the test can return is a third -- no matter how large the
# observed effect is. Reporting the count of attainable arrangements alone hides
# that, so the resolution is computed here alongside the p-value.
#
# `minimum_attainable_p` is the p-value the test would return if the observed
# statistic were the largest the permutation distribution can produce.
#' @keywords internal
#' @noRd
.distinct_values <- function(x, tol) {
  x <- sort(x[is.finite(x)])
  if (!length(x)) return(0L)
  1L + sum(diff(x) > tol * max(1, max(abs(x))))
}

# `alpha` is the level the resolution note is written against: a design whose
# smallest attainable p-value exceeds it cannot reject at that level however
# large the effect. It is a reporting threshold only -- `distinct_statistics`
# and `minimum_attainable_p` are recorded regardless -- so it is taken from the
# `distmer.resolution_alpha` option rather than threaded through every engine.
#' @keywords internal
#' @noRd
.permutation_inference <- function(statistics, observed, exact, tol = 1e-10,
                                   alpha = getOption("distmer.resolution_alpha",
                                                     0.05)) {
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`distmer.resolution_alpha` must be one number strictly between ",
         "zero and one.", call. = FALSE)
  }
  count <- length(statistics)
  empty <- list(
    p_value = NA_real_, distinct_statistics = 0L,
    minimum_attainable_p = NA_real_, resolution_note = character()
  )
  if (!count) return(empty)

  # Deciding whether a permuted statistic ties with the observed one is a
  # numerical-equality question, so the slack has to cover both the magnitude of
  # the statistic and the error accumulated in computing it. A pseudo-F is a
  # ratio of traces of projected Gower matrices; under strong separation that
  # chain loses several digits, and gaps of ~1e-9 relative are ordinary. An
  # absolute 1e-10 tolerance then fails to count the identity permutation as
  # extreme and an exact test reports p = 0, which no permutation test can
  # legitimately produce. `sqrt(.Machine$double.eps)` is the floor used for
  # numerical-equality checks elsewhere in the package; erring loose here is the
  # conservative direction, since it can only keep genuine ties tied.
  slack <- max(100 * tol, sqrt(.Machine$double.eps)) * max(1, abs(observed))
  extreme <- sum(statistics >= observed - slack, na.rm = TRUE)
  # Under exact enumeration the identity permutation is a member of the set and
  # reproduces the observed statistic, so at least one arrangement always ties
  # and p >= 1/N holds as a theorem. The observed statistic and the permuted
  # ones are reached by different arithmetic paths, though, and under severe
  # conditioning (a pseudo-F of 1e12, say) that tie can be lost to rounding by
  # more than any fixed relative slack. Enforcing the bound restores the
  # guarantee without depending on how well-conditioned the design is.
  if (isTRUE(exact)) extreme <- max(extreme, 1L)
  p_value <- if (isTRUE(exact)) extreme / count else (extreme + 1) / (count + 1)

  finite <- statistics[is.finite(statistics)]
  if (!length(finite)) {
    empty$p_value <- p_value
    return(empty)
  }
  distinct <- .distinct_values(finite, tol)
  ties <- sum(finite >= max(finite) -
                max(100 * tol, sqrt(.Machine$double.eps)) *
                  max(1, abs(max(finite))))
  minimum <- if (isTRUE(exact)) ties / count else 1 / (count + 1)

  note <- character()
  if (is.finite(minimum) && minimum > alpha) {
    note <- paste0(
      "The permutation distribution holds only ", distinct,
      " distinct statistic(s); the smallest p-value this design can return is ",
      format(round(minimum, 4)), ", so it cannot reach ", alpha,
      " however large the effect."
    )
  }
  list(p_value = p_value, distinct_statistics = distinct,
       minimum_attainable_p = minimum, resolution_note = note)
}

# Is a mean square numerically zero?
#
# A mean square carries the squared units of the dissimilarity, so testing it
# against an absolute floor makes the answer depend on how the response happens
# to be scaled: the same analysis passes on a response of ordinary magnitude and
# fails after multiplying it by 1e-4. Total inertia sets the scale instead, so
# the question becomes "is this stratum negligible relative to the variation in
# the data", which is what the guard is actually for and is scale free.
#' @keywords internal
#' @noRd
.mean_square_is_zero <- function(ms, G = NULL, tol = 1e-10, reference = NULL) {
  if (length(ms) != 1L || !is.finite(ms)) return(TRUE)
  # Total *absolute* inertia, not the signed trace: a retained non-Euclidean
  # spectrum can carry negative eigenvalues large enough to make the trace
  # negative, and that is a supported result rather than an absence of
  # variation. Magnitudes scale with the response exactly as the mean square
  # does, so the ratio stays scale free either way.
  if (is.null(reference)) {
    reference <- sum(abs(diag(G))) / max(1L, nrow(G) - 1L)
  }
  reference <- abs(reference)
  # No variation anywhere leaves nothing to test against.
  if (!is.finite(reference) || reference <= 0) return(TRUE)
  abs(ms) <= tol * reference
}
