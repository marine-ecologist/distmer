# PRIMER-style pseudo-F construction from expected mean squares.
#
# `quasi_denominator()` finds a linear combination of candidate mean squares
# whose EMS matches the numerator's. Under an unbalanced design that solution
# can carry a negative coefficient on the residual, which is not a usable mean
# square. PRIMER avoids this by balancing the residual across both sides of the
# ratio instead of forcing it into the denominator:
#
#   numerator EMS under H0 = 1*V(Res) + sum_k a_k V_k
#   stratum k EMS          = 1*V(Res) + A_k V_k
#
# Matching each random component exactly needs c_k = a_k / A_k, and the strata
# then carry s = sum_k c_k units of V(Res). When s exceeds one the numerator is
# short of residual and takes (s - 1) * MS_Res alongside the term; when s falls
# short the denominator takes (1 - s) * MS_Res. Either way every coefficient
# stays non-negative, and Satterthwaite degrees of freedom apply to whichever
# side became a combination -- so a numerator df need not be a whole number.
#
# Verified against PRIMER's own printed construction for an unbalanced
# three-stratum reef design: pseudo-F, numerator df and denominator df all
# reproduce to every published digit.
#' @keywords internal
#' @noRd
.satterthwaite <- function(coefficients, mean_squares, degrees) {
  keep <- abs(coefficients) > 0
  contribution <- coefficients[keep] * mean_squares[keep]
  total <- sum(contribution)
  if (!is.finite(total) || total <= 0) {
    return(list(value = total, df = NA_real_))
  }
  list(
    value = total,
    df = total^2 / sum(contribution^2 / degrees[keep])
  )
}

#' @keywords internal
#' @noRd
.ems_f_construction <- function(numerator_ems, numerator_ms, numerator_df,
                                candidate_ems, candidate_ms, candidate_df,
                                residual = "residual", after = NULL,
                                tol = 1e-10) {
  components <- setdiff(names(numerator_ems), residual)
  # The null leaves nothing but the residual component when every other
  # coefficient is zero -- either because no other component exists, or because
  # the term's own component was the only one and testing it sets it to zero.
  exhausted <- !length(components) ||
    all(abs(numerator_ems[components]) <= tol)
  if (exhausted) {
    # Nothing but the residual survives, so the denominator carries the residual
    # component alone. That does not make it the residual stratum: in a nested
    # design an intermediate stratum can have exactly the same expected mean
    # square, and the test belongs against the stratum immediately below the
    # term rather than against the bottom of the hierarchy. Sponge within
    # chamber within treatment and the residual can both read 1*V(Res), and
    # chamber is tested over sponge, not over the residual -- different degrees
    # of freedom and a different answer.
    #
    # `after` names the term's own position among the candidates, which are
    # ordered coarsest to finest; the denominator is then the first stratum
    # below it that carries no other component. Without it the residual is the
    # only defensible choice.
    row <- integer()
    if (!is.null(after)) {
      position <- match(after, rownames(candidate_ems))
      if (!is.na(position) && position < nrow(candidate_ems)) {
        below <- seq.int(position + 1L, nrow(candidate_ems))
        empty <- below[vapply(below, function(k) {
          all(abs(candidate_ems[k, , drop = TRUE][
            setdiff(colnames(candidate_ems), residual)]) <= tol)
        }, logical(1))]
        if (length(empty)) row <- empty[[1L]]
      }
    }
    if (!length(row)) row <- which(rownames(candidate_ems) == "Residual")
    if (!length(row)) row <- nrow(candidate_ems)
    row <- row[[length(row)]]
    coefficient <- numeric(nrow(candidate_ems))
    coefficient[[row]] <- 1
    return(list(
      statistic = numerator_ms / candidate_ms[[row]],
      numerator_ms = numerator_ms, numerator_df = numerator_df,
      denominator_ms = candidate_ms[[row]], denominator_df = candidate_df[[row]],
      coefficients = coefficient, numerator_residual = 0,
      strata_used = rownames(candidate_ems)[row], single = TRUE
    ))
  }
  # The residual stratum matches nothing but the residual, so it is the balance
  # term rather than a candidate for the component match.
  # More than one stratum can carry the residual component alone: in a nested
  # design an intermediate stratum may have exactly the same expected mean
  # square as the residual. Candidates run coarsest to finest, so the balance
  # term is the last of them; the earlier ones remain available as denominators
  # for the terms above them.
  residual_row <- which(
    apply(candidate_ems[, components, drop = FALSE], 1L,
          function(row) all(abs(row) <= tol))
  )
  if (!length(residual_row)) {
    return(NULL)
  }
  residual_row <- residual_row[[length(residual_row)]]
  strata <- setdiff(seq_len(nrow(candidate_ems)), residual_row)
  # Solve for the stratum weights that reproduce the numerator's random
  # components exactly. Orthogonalising the nested components makes the system
  # triangular, but a stratum can still carry a small amount of a component
  # below it under imbalance, so this is a solve rather than a lookup.
  # A stratum carrying no component at all contributes a zero column and would
  # make the system singular; it can never be part of a match, so leave it out.
  informative <- vapply(strata, function(k) {
    any(abs(candidate_ems[k, components, drop = TRUE]) > tol)
  }, logical(1))
  strata <- strata[informative]
  if (!length(strata)) return(NULL)
  A <- t(candidate_ems[strata, components, drop = FALSE])
  target <- numerator_ems[components]
  weights <- tryCatch(
    if (nrow(A) == ncol(A)) solve(A, target) else .pseudo_inverse_solve(A, target, tol = tol),
    error = function(e) NULL
  )
  if (is.null(weights) || any(!is.finite(weights))) return(NULL)
  if (max(abs(A %*% weights - target)) > sqrt(tol) * max(1, max(abs(target)))) {
    return(NULL)
  }
  weights[abs(weights) <= tol] <- 0
  if (any(weights < -tol)) return(NULL)
  coefficient <- numeric(nrow(candidate_ems))
  coefficient[strata] <- weights
  carried <- sum(coefficient[strata])

  if (carried > 1 + tol) {
    numerator <- .satterthwaite(
      c(carried - 1, 1), c(candidate_ms[[residual_row]], numerator_ms),
      c(candidate_df[[residual_row]], numerator_df))
    denominator <- .satterthwaite(coefficient[strata], candidate_ms[strata],
                                  candidate_df[strata])
    numerator_extra <- carried - 1
  } else {
    coefficient[[residual_row]] <- 1 - carried
    numerator <- list(value = numerator_ms, df = numerator_df)
    denominator <- .satterthwaite(coefficient, candidate_ms, candidate_df)
    numerator_extra <- 0
  }
  if (!is.finite(denominator$value) || denominator$value <= 0) return(NULL)
  used <- rownames(candidate_ems)[abs(coefficient) > tol]
  list(
    statistic = numerator$value / denominator$value,
    numerator_ms = numerator$value, numerator_df = numerator$df,
    denominator_ms = denominator$value, denominator_df = denominator$df,
    coefficients = coefficient, numerator_residual = numerator_extra,
    strata_used = used,
    single = length(used) == 1L
  )
}

# Restricted-model expected mean squares.
#
# PRIMER states its convention in the workspace header: "Fixed effects sum to
# zero for mixed terms". Under that restricted model the mean square of a random
# stratum carries no variance component from a random-by-fixed interaction
# nested inside it -- E(MS_ISLAND) contains no sigma^2_(ISLAND x DEPTH), even
# though E(MS_(ISLAND x DEPTH)) does. Only the stratum rows are restricted; the
# contrast's own expected mean square keeps every component it carries.
#
# The zero belongs in the accounting, not in the component. Getting it by
# orthogonalising Z against the containing stratum instead would change that
# component's covariance everywhere it appears -- including in the contrast's
# own EMS -- and leave the inner coefficient too small.
.restricted_ems <- function(candidate, interaction_terms) {
  components <- colnames(candidate)
  parts <- strsplit(components, ":", fixed = TRUE)
  names(parts) <- components
  for (inner in intersect(interaction_terms, components)) {
    for (stratum in rownames(candidate)) {
      if (identical(stratum, inner) || !stratum %in% components) next
      if (all(parts[[stratum]] %in% parts[[inner]])) {
        candidate[stratum, inner] <- 0
      }
    }
  }
  candidate
}
