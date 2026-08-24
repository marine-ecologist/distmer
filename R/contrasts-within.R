# Pairwise contrasts of one factor within the cells of others.
#
# The denominator for such a contrast cannot be inherited from the parent fit.
# Inside a single cell the lower strata often collapse -- an ecoregion sampled
# on two islands leaves the unit-by-sub-plot stratum with no degrees of freedom
# at all -- so the parent's denominator does not exist there. PRIMER derives a
# denominator per pair by expected-mean-square matching against whichever strata
# survive inside that subset, which is why three pairs within one cell can carry
# three different denominators and three different Satterthwaite df.
#
# This reproduces that: subset, rebuild the declared random strata on the
# subset, drop the ones with no degrees of freedom, and match the contrast's EMS
# against the survivors with `.ems_f_construction()`.

# Sequentially-built strata for a subset, given the declared grouping factors
# (outermost first) and the fixed columns the contrast itself spans.
.within_strata <- function(G, groups, fixed_columns, tol = 1e-10) {
  n <- nrow(G)
  # Every stratum is built with the contrast already in the model, exactly as
  # the parent engine builds its sub-plot stratum after the sub-plot fixed
  # columns. Leaving the contrast out lets the deepest grouping absorb it, which
  # inflates that stratum and empties the residual.
  base <- cbind(matrix(1, n, 1L), fixed_columns)
  previous <- projection_matrix(base, tol = tol)
  spanned <- base
  strata <- list()
  for (name in names(groups)) {
    level <- droplevels(factor(groups[[name]]))
    if (nlevels(level) < 2L) next
    Z <- stats::model.matrix(~ 0 + level)
    current <- projection_matrix(cbind(spanned, Z), tol = tol)
    df <- current$rank - previous$rank
    if (df > 0L) {
      strata[[name]] <- list(H = current$H - previous$H, df = df,
                             ms = trace_hgh(current$H - previous$H, G) / df,
                             Z = Z)
    }
    previous <- current
    spanned <- cbind(spanned, Z)
  }
  # The residual is whatever the deepest grouping leaves.
  full <- projection_matrix(spanned, tol = tol)
  residual_df <- n - full$rank
  if (residual_df > 0L) {
    H <- diag(n) - full$H
    strata[["Residual"]] <- list(H = H, df = residual_df,
                                 ms = trace_hgh(H, G) / residual_df, Z = NULL)
  }
  strata
}

# One pairwise contrast inside one cell.
.within_pair <- function(G, keep, values, pair, groups, permutations, seed,
                         factor_name = NULL, tol = 1e-10) {
  index <- which(keep & values %in% pair)
  if (length(index) < 3L) return(NULL)
  Gs <- G[index, index, drop = FALSE]
  Gs <- gower_center(stats::as.dist(sqrt(pmax(
    outer(diag(Gs), diag(Gs), "+") - 2 * Gs, 0))))
  level <- droplevels(factor(values[index]))
  if (nlevels(level) != 2L) return(NULL)
  n <- length(index)
  sub_groups <- lapply(groups, function(g) droplevels(factor(g[index])))
  names(sub_groups) <- names(groups)

  contrast <- ifelse(level == pair[[1L]],
                     1 / sum(level == pair[[1L]]),
                     -1 / sum(level == pair[[2L]]))
  H_num <- tcrossprod(contrast) / sum(contrast^2)
  numerator_ms <- trace_hgh(H_num, Gs)

  strata <- .within_strata(Gs, sub_groups, stats::model.matrix(~ level)[, -1L, drop = FALSE],
                           tol = tol)
  if (!length(strata)) return(NULL)

  # EMS of the contrast and of every surviving stratum, on the subset's own
  # covariance bases.
  bases <- list(residual = diag(n))
  for (name in names(strata)) {
    if (!is.null(strata[[name]]$Z)) bases[[name]] <- tcrossprod(strata[[name]]$Z)
  }
  # Each component's covariance is Z Z' as declared -- a nested component is NOT
  # orthogonalised against the one containing it. Orthogonalising would force
  # every component to be carried by exactly one stratum, which is not how
  # expected mean squares behave: an inner component legitimately appears in the
  # outer stratum's EMS as well, and the matching in .ems_f_construction() is
  # what resolves it. Orthogonalising here made the inner coefficient too small
  # by up to 18%, worst where replication was most unbalanced.
  ems_of <- function(H, df) {
    vapply(bases, function(V) trace_product(H, V) / df, numeric(1))
  }
  candidate <- do.call(rbind, lapply(names(strata), function(name) {
    ems_of(strata[[name]]$H, strata[[name]]$df)
  }))
  rownames(candidate) <- names(strata)
  # Strata whose grouping involves the contrasted factor are random-by-fixed
  # interactions, so the restricted model zeroes them out of any stratum they
  # sit inside.
  interaction_terms <- if (is.null(factor_name)) character() else {
    Filter(function(name) factor_name %in% strsplit(name, ":", fixed = TRUE)[[1L]],
           names(strata))
  }
  candidate <- .restricted_ems(candidate, interaction_terms)
  candidate_ms <- vapply(strata, `[[`, numeric(1), "ms")
  candidate_df <- vapply(strata, `[[`, numeric(1), "df")

  numerator_ems <- ems_of(H_num, 1L)
  numerator_ems[["residual"]] <- 1
  built <- .ems_f_construction(numerator_ems, numerator_ms, 1L,
                               candidate, candidate_ms, candidate_df, tol = tol)
  if (is.null(built)) return(NULL)

  list(index = index, statistic = built$statistic,
       numerator_residual = built$numerator_residual,
       residual_ms = if ("Residual" %in% names(strata)) {
         strata[["Residual"]]$ms
       } else 0,
       t = sqrt(max(built$statistic, 0)),
       numerator_ms = numerator_ms,
       denominator_ms = built$denominator_ms,
       denominator_df = built$denominator_df,
       # A weight that survives the solve at machine epsilon is zero: it should
       # not be named as a stratum the test leans on.
       denominator = paste(sprintf("%.5g*%s",
                                   built$coefficients[built$coefficients > tol],
                                   rownames(candidate)[built$coefficients > tol]),
                           collapse = " + "),
       strata = strata, H_num = H_num, G = Gs, groups = sub_groups,
       coefficients = built$coefficients)
}

# Permutation at the level of a grouping stratum.
#
# When the denominator's finest term is a grouping stratum rather than the
# residual, the exchangeable entity is a whole unit of that stratum -- an island
# x depth cell here -- not an observation. Permuting observations inside the
# island instead treats each of the ~40 samples in a cell as an independent
# replicate, which is pseudo-replication: the null comes out far too tight and
# the p-value too small. PRIMER's "Unique perms" column shows the intended
# space directly -- every cell with six occupied units reports exactly
# 720 = 6! -- so units are permuted freely, and because units carry different
# numbers of samples all u! orderings give distinct statistics.
#
# Permuting whole units of unequal size is not a permutation of the rows, so it
# is done on the unit centroids. That costs nothing: writing L for the operator
# that replaces every observation by its unit centroid, any projector whose
# column space is constant within units -- the contrast, and every stratum
# coarser than or equal to the unit -- satisfies
#
#   tr(H L G L') = tr(H G),
#
# so the observed statistic is identical either way and only the null changes.
# Only the residual stratum is altered by L, and it carries no weight whenever
# this path is taken.
.within_unit_permutation <- function(pair, finest, positive, permutations,
                                     max_exact, tol = 1e-10) {
  unit <- droplevels(factor(pair$groups[[finest]]))
  A <- stats::model.matrix(~ 0 + unit)
  sizes <- colSums(A)
  units <- ncol(A)
  # Gower matrix of the unit centroids.
  inverse <- 1 / sizes
  G_unit <- crossprod(A, pair$G %*% A) * outer(inverse, inverse)
  kernel <- function(H) crossprod(A, H %*% A)
  K_numerator <- kernel(pair$H_num)
  K_strata <- lapply(positive, function(k) kernel(pair$strata[[k]]$H))
  weights <- pair$coefficients[positive]
  divisors <- vapply(pair$strata[positive], `[[`, numeric(1), "df")

  # Freedman-Lane, as PRIMER states it does: residuals are permuted under the
  # reduced model, which holds every stratum coarser than the unit. Because
  # units carry different numbers of samples the reduced-model projector is
  # weighted by unit size -- it is the observation-level fit expressed on unit
  # centroids, not an unweighted average of them.
  coarser <- names(pair$strata)[seq_len(match(finest, names(pair$strata)) - 1L)]
  first_of <- function(values) {
    droplevels(factor(vapply(split(as.character(values), unit),
                             `[[`, character(1), 1L)))
  }
  reduced_X <- matrix(1, units, 1L)
  for (name in coarser) {
    reduced_X <- cbind(reduced_X,
                       stats::model.matrix(~ 0 + first_of(pair$groups[[name]])))
  }
  # Weighted least squares through a symmetric eigendecomposition, so the
  # rank-deficient indicator coding (intercept plus a full set of dummies) is
  # handled without dropping columns by hand.
  weighted <- crossprod(reduced_X, reduced_X * sizes)
  spectrum <- eigen(weighted, symmetric = TRUE)
  keep <- spectrum$values > max(tol, spectrum$values[[1L]] * tol)
  inverse <- spectrum$vectors[, keep, drop = FALSE] %*%
    (t(spectrum$vectors[, keep, drop = FALSE]) / spectrum$values[keep])
  P_reduced <- reduced_X %*% (inverse %*% t(reduced_X * sizes))
  R_reduced <- diag(units) - P_reduced

  plan <- .permutation_orders(units, permutations, max_exact)
  if (!nrow(plan$orders)) {
    return(list(p_value = NA_real_, used = 0L, exact = FALSE,
                distinct_statistics = 0L, minimum_attainable_p = NA_real_,
                resolution_note = character()))
  }
  # When the matched weights carry more than one unit of residual,
  # .ems_f_construction() moves the surplus onto the numerator rather than
  # allowing a negative weight, so the statistic is
  #
  #   (numerator MS + surplus * residual MS) / denominator.
  #
  # The surplus term is constant here: whole units move intact, so the residual
  # sum of squares within units -- and its degrees of freedom -- do not change.
  # Dropping it would leave the null on a different scale from the observed
  # statistic it is compared against.
  surplus <- pair$numerator_residual * pair$residual_ms
  statistic <- function(index) {
    shifted <- P_reduced + R_reduced[index, , drop = FALSE]
    permuted <- shifted %*% G_unit %*% t(shifted)
    denominator <- sum(vapply(seq_along(K_strata), function(k) {
      weights[[k]] * sum(K_strata[[k]] * permuted) / divisors[[k]]
    }, numeric(1)))
    (sum(K_numerator * permuted) + surplus) / denominator
  }
  # The observed configuration must reproduce the statistic the table reports.
  identity_statistic <- statistic(seq_len(units))
  if (!isTRUE(all.equal(identity_statistic, pair$statistic,
                        tolerance = 1e-6))) {
    stop("Unit-level permutation does not reproduce the observed statistic ",
         "(", format(identity_statistic), " vs ", format(pair$statistic), ").",
         call. = FALSE)
  }
  statistics <- .permuted_statistics(nrow(plan$orders),
                                     function(row) statistic(plan$orders[row, ]))
  inference <- .permutation_inference(statistics, pair$statistic, plan$exact, tol)
  c(inference, list(used = nrow(plan$orders),
                    exact = plan$exact && nrow(plan$orders) > 0L))
}

# Permutation distribution for one within-cell contrast. The reduced model holds
# whatever grouping survives in the subset, and residuals move only inside those
# groups; with no grouping left the exchange is unrestricted. The denominator is
# recombined under each permutation with the coefficients used for the observed
# statistic, so the null matches the test.
.within_permutation <- function(pair, permutations, max_exact, tol = 1e-10) {
  n <- nrow(pair$G)
  # The exchangeable entity is a unit of whichever stratum sits deepest in the
  # denominator. `strata` runs outermost to innermost with the residual last, so
  # the last stratum carrying weight names it.
  positive <- which(pair$coefficients > tol)
  if (length(positive)) {
    finest <- names(pair$strata)[[max(positive)]]
    if (!identical(finest, "Residual")) {
      return(.within_unit_permutation(pair, finest, positive, permutations,
                                      max_exact, tol = tol))
    }
  }
  # The residual carries the denominator, so observations are the units. They
  # still move only inside the OUTERMOST surviving grouping, so an observation
  # crosses the levels being contrasted without leaving its experimental unit.
  grouping <- Filter(function(name) !identical(name, "Residual"), names(pair$strata))
  deepest <- if (length(grouping)) pair$groups[[grouping[[1L]]]] else NULL
  reduced_H <- if (is.null(deepest)) {
    projection_matrix(matrix(1, n, 1L), tol = tol)$H
  } else {
    projection_matrix(stats::model.matrix(~ 0 + deepest), tol = tol)$H
  }
  plan <- if (is.null(deepest)) {
    .permutation_orders(n, permutations, max_exact)
  } else {
    .stratified_permutation_plan(deepest, permutations, max_exact)
  }
  if (!nrow(plan$orders)) {
    return(list(p_value = NA_real_, used = 0L, exact = FALSE,
                distinct_statistics = 0L, minimum_attainable_p = NA_real_,
                resolution_note = character()))
  }
  positive <- pair$coefficients > tol
  statistics <- if (sum(positive) == 1L && pair$numerator_residual == 0) {
    .freedman_lane_statistics(
      pair$G, reduced_H, pair$H_num, 1L,
      pair$strata[[which(positive)]]$H, pair$strata[[which(positive)]]$df,
      plan$orders, tol = tol)
  } else {
    .freedman_lane_quasi_statistics(
      pair$G, reduced_H, pair$H_num, 1L,
      lapply(pair$strata, `[[`, "H"),
      vapply(pair$strata, `[[`, numeric(1), "df"),
      pair$coefficients, plan$orders,
      numerator_residual = pair$numerator_residual,
      residual_index = match("Residual", names(pair$strata)), tol = tol)
  }
  inference <- .permutation_inference(statistics, pair$statistic, plan$exact, tol)
  c(inference, list(used = nrow(plan$orders),
                    exact = plan$exact && nrow(plan$orders) > 0L))
}

#' Pairwise contrasts within the cells of other factors
#'
#' Contrasts the levels of one fixed factor separately inside every combination
#' of the `within` factors, as PRIMER's "within level" pairwise tables do. Each
#' pair gets its own denominator: the declared random strata are rebuilt on that
#' subset, any left with no degrees of freedom are dropped, and the contrast's
#' expected mean square is matched against whichever survive. Three pairs inside
#' one cell can therefore carry three different denominators and three different
#' Satterthwaite degrees of freedom.
#'
#' @param fit A `permanova` object.
#' @param factor Name of the fixed factor whose levels are contrasted.
#' @param within Character vector of factors whose combinations define the cells.
#' @param permutations Permutations per contrast.
#' @param seed Optional random seed.
#' @param adjust Multiplicity method passed to [stats::p.adjust()], applied
#'   within each cell rather than across the whole table.
#' @param max_exact Maximum permutation space allowed for exact enumeration.
#' @param tol Numerical tolerance.
#' @return A data frame of contrasts, one row per pair per cell.
#' @export
permanova_contrasts_within <- function(fit, factor, within,
                                       permutations = 999L, seed = NULL,
                                       adjust = "none", max_exact = 100000L,
                                       tol = 1e-10) {
  if (!inherits(fit, "permanova")) .stopf("`fit` must be a permanova() object")
  design <- fit$design %||% .stopf("fit has no $design")
  data <- design$data
  needed <- c(factor, within)
  missing <- needed[!needed %in% names(data)]
  if (length(missing)) {
    .stopf("not in fit$design$data: %s", paste(missing, collapse = ", "))
  }
  if (length(factor) != 1L) .stopf("`factor` must name one fixed factor")
  if (!is.null(seed)) set.seed(seed)

  values <- droplevels(factor(data[[factor]]))
  groups <- lapply(design$random_terms, function(term) {
    droplevels(factor(design$groups[[term$label]]))
  })
  names(groups) <- vapply(design$random_terms, `[[`, character(1), "label")
  G <- design$G %||% .stopf("fit$design has no $G")

  cells <- unique(data[, within, drop = FALSE])
  cells <- cells[do.call(order, as.list(cells)), , drop = FALSE]
  rows <- list()
  for (i in seq_len(nrow(cells))) {
    keep <- rep(TRUE, nrow(data))
    for (v in within) keep <- keep & data[[v]] == cells[[v]][i]
    present <- levels(droplevels(values[keep]))
    if (length(present) < 2L) next
    pairs <- utils::combn(present, 2L, simplify = FALSE)
    cell_rows <- list()
    for (pr in pairs) {
      built <- .within_pair(G, keep, values, pr, groups, permutations, seed,
                            factor_name = factor, tol = tol)
      if (is.null(built)) next
      perm <- .within_permutation(built, permutations, max_exact, tol = tol)
      cell_rows[[length(cell_rows) + 1L]] <- data.frame(
        cells[i, , drop = FALSE],
        contrast = paste(pr, collapse = " - "),
        t = built$t, pseudo_F = built$statistic,
        denominator = built$denominator,
        denominator_df = built$denominator_df,
        p_perm = perm$p_value, unique_perms = perm$used,
        distinct_statistics = perm$distinct_statistics,
        minimum_attainable_p = perm$minimum_attainable_p,
        row.names = NULL, stringsAsFactors = FALSE
      )
    }
    if (!length(cell_rows)) next
    block <- do.call(rbind, cell_rows)
    block$p_adjusted <- stats::p.adjust(block$p_perm, method = adjust)
    rows[[length(rows) + 1L]] <- block
  }
  if (!length(rows)) .stopf("no contrast was estimable in any cell")
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  attr(out, "within") <- within
  attr(out, "factor") <- factor
  class(out) <- c("permanova_contrasts_within", "data.frame")
  out
}

#' @param x A `permanova_contrasts_within` object.
#' @param digits Significant digits for the statistic columns.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @noRd
#' @exportS3Method base::print
print.permanova_contrasts_within <- function(x, digits = 4, ...) {
  cells <- attr(x, "within")
  contrasted <- attr(x, "factor")
  frame <- as.data.frame(x)
  # A subset may have dropped the grouping columns; fall back to a flat table
  # rather than failing to print.
  if (is.null(cells) || !all(cells %in% names(frame))) {
    return(print(frame, ...))
  }
  cat(sprintf("Pairwise contrasts of '%s' within %s\n",
              contrasted %||% "factor", paste(cells, collapse = " x ")))
  cat(sprintf("  %d contrasts in %d cells; strata rebuilt per pair\n\n",
              nrow(frame), nrow(unique(frame[, cells, drop = FALSE]))))
  keys <- do.call(paste, c(frame[, cells, drop = FALSE], list(sep = " | ")))
  for (key in unique(keys)) {
    block <- frame[keys == key, , drop = FALSE]
    cat(paste0(key, "\n"))
    show <- data.frame(
      contrast = block$contrast,
      t = round(block$t, digits),
      df = round(block$denominator_df, 2),
      denominator = block$denominator,
      p_perm = ifelse(is.na(block$p_perm), "-",
                      format.pval(block$p_perm, digits = 3, eps = 1e-4)),
      check.names = FALSE, stringsAsFactors = FALSE
    )
    print(show, row.names = FALSE, right = FALSE)
    cat("\n")
  }
  invisible(x)
}
