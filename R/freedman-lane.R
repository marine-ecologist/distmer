# Shared Freedman-Lane permutation kernel and permutation progress reporting.
#
# Every Freedman-Lane engine in this package forms
#
#   A = H_reduced + Pi %*% R,      R = I - H_reduced,
#
# for a permutation matrix Pi, and then needs tr(H A G A') for the numerator
# and denominator projectors. Forming `A %*% G %*% t(A)` explicitly costs two
# dense n-by-n-by-n products for every permutation. Expanding the product,
#
#   A G A' = Hr G Hr + Hr G R Pi' + Pi R G Hr + Pi (R G R) Pi',
#
# each piece collapses against a fixed matrix:
#
#   tr(H Hr G Hr)      is constant across permutations;
#   tr(H Hr G R Pi')   and tr(H Pi R G Hr) are equal for symmetric H and read
#                      n entries out of the fixed matrix `base = R G Hr`;
#   tr(H Pi M Pi')     is the only whole-matrix term, with M = R G R fixed.
#
# The remaining term is where the low-rank structure is spent. A hypothesis
# projector carries rank equal to its degrees of freedom -- typically a handful
# -- so writing H = W diag(s) W' with W an n-by-r basis gives
#
#   tr(H Pi M Pi') = tr(s W' Pi M Pi' W) = sum_k s_k * V[, k]' M V[, k],
#
# for V = Pi' W, a cheap row re-indexing of W. Each permutation therefore costs
# one `M %*% V` product with r columns instead of two dense n-by-n products.
# Residual projectors have rank close to n, so those are represented through
# their complement, I - H, using tr(Pi M Pi') = tr(M): permutation similarity
# leaves the total unchanged, and only the low-rank correction is recomputed.
#
# Set-up costs one symmetric eigendecomposition per projector plus a fixed
# handful of dense products, repaid within a few permutations.

# Represent a symmetric matrix as `scale * I + W diag(signs) W'` with W as
# narrow as possible: directly when the rank is low, through the complement
# when it is high. Returns NULL when neither representation is narrow, leaving
# the caller on the dense path.
.low_rank_representation <- function(H, tol = 1e-10) {
  n <- nrow(H)
  if (max(abs(H - t(H))) > tol) {
    return(NULL)
  }
  spectrum <- eigen(H, symmetric = TRUE)
  values <- spectrum$values
  width <- n %/% 2L
  direct <- abs(values) > tol
  complement <- abs(values - 1) > tol
  scale <- if (sum(direct) <= sum(complement)) 0 else 1
  keep <- if (scale == 0) direct else complement
  if (sum(keep) > width) {
    return(NULL)
  }
  retained <- values[keep] - scale
  basis <- spectrum$vectors[, keep, drop = FALSE]
  list(
    scale = scale,
    signs = sign(retained),
    W = sweep(basis, 2L, sqrt(abs(retained)), `*`)
  )
}

.freedman_lane_kernel <- function(G, H_reduced, projectors, tol = 1e-10) {
  n <- nrow(G)
  residual <- diag(n) - H_reduced
  M <- residual %*% (G %*% residual)
  M <- (M + t(M)) / 2
  base <- residual %*% (G %*% H_reduced)
  reduced_G_reduced <- H_reduced %*% (G %*% H_reduced)
  trace_M <- sum(diag(M))
  seq_n <- seq_len(n)
  row_offset <- (seq_n - 1L) * n

  parts <- lapply(projectors, function(H) {
    part <- list(
      H = H,
      constant = sum(H * t(reduced_G_reduced)),
      rank = .low_rank_representation(H, tol = tol)
    )
    if (is.null(part$rank)) {
      # Dense fallback: still O(n^2) per permutation, just not in BLAS.
      part$right <- base %*% H
      part$left <- if (max(abs(H - t(H))) <= tol) NULL else H %*% t(base)
    } else {
      part$base_W <- base %*% part$rank$W
      part$signed_W <- sweep(part$rank$W, 2L, part$rank$signs, `*`)
    }
    part
  })

  # One `M %*% V` covering every low-rank projector at once, so a permutation
  # makes a single pass over M no matter how many hypotheses share it.
  widths <- vapply(parts, function(part) {
    if (is.null(part$rank)) 0L else ncol(part$rank$W)
  }, integer(1))
  starts <- cumsum(c(0L, widths))
  blocks <- lapply(seq_along(widths), function(position) {
    seq_len(widths[[position]]) + starts[[position]]
  })
  stacked <- if (sum(widths)) {
    do.call(cbind, lapply(parts[widths > 0L], function(part) part$rank$W))
  } else NULL

  evaluate <- function(index) {
    values <- numeric(length(parts))
    quadratic <- numeric(0)
    if (!is.null(stacked)) {
      inverse <- integer(n)
      inverse[index] <- seq_n
      V <- stacked[inverse, , drop = FALSE]
      quadratic <- colSums(V * (M %*% V))
    }
    for (position in seq_along(parts)) {
      part <- parts[[position]]
      if (is.null(part$rank)) {
        permuted <- M[index, index, drop = FALSE]
        cross <- sum(part$right[index + row_offset])
        cross <- if (is.null(part$left)) {
          2 * cross
        } else {
          cross + sum(part$left[seq_n + (index - 1L) * n])
        }
        # `permuted` is symmetric, so sum(H * permuted) is tr(H permuted).
        values[[position]] <- part$constant + cross + sum(part$H * permuted)
        next
      }
      whole <- part$rank$scale
      cross <- sum(part$base_W[index, , drop = FALSE] * part$signed_W)
      if (whole != 0) {
        cross <- cross + whole * sum(base[index + row_offset])
      }
      values[[position]] <- part$constant + 2 * cross + whole * trace_M +
        sum(part$rank$signs * quadratic[blocks[[position]]])
    }
    values
  }

  list(n = n, evaluate = evaluate)
}

# Drive a complete Freedman-Lane permutation distribution for one pseudo-F.
#
# `orders` holds one permutation per row. `index` optionally maps a row of
# `orders` onto an observation-level index vector, which is how the engines
# that permute intact blocks expand a unit permutation.
.freedman_lane_statistics <- function(G, H_reduced, H_numerator, numerator_df,
                                      H_denominator, denominator_df, orders,
                                      index = NULL, label = NULL, tol = 1e-10) {
  total <- nrow(orders)
  statistics <- numeric(total)
  if (!total) {
    return(statistics)
  }
  kernel <- .freedman_lane_kernel(
    G, H_reduced, list(H_numerator, H_denominator), tol = tol
  )
  progress <- .permutation_progress(total, label)
  on.exit(progress$close(), add = TRUE)
  for (row in seq_len(total)) {
    permutation <- orders[row, ]
    if (!is.null(index)) {
      permutation <- index(permutation)
    }
    traces <- kernel$evaluate(permutation)
    statistics[[row]] <- (traces[[1L]] / numerator_df) /
      (traces[[2L]] / denominator_df)
    progress$tick(row)
  }
  statistics
}

#' Report permutation progress
#'
#' Returns a `tick`/`close` pair driving a text progress bar. Reporting is off
#' unless the session is interactive, and is suppressed for short runs where a
#' bar is noise rather than information. Set `options(distmer.progress = TRUE)`
#' to force it on (for example inside a script), `FALSE` to silence it, and
#' `distmer.progress_min` to change the length at which a run is worth
#' reporting.
#'
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @keywords internal
#' @noRd
.permutation_progress <- function(total, label = NULL) {
  quiet <- list(
    tick = function(iteration) invisible(NULL),
    close = function() invisible(NULL)
  )
  enabled <- isTRUE(getOption("distmer.progress", interactive())) &&
    length(total) == 1L && is.finite(total) &&
    total >= getOption("distmer.progress_min", 200L)
  if (!enabled) {
    return(quiet)
  }
  if (!is.null(label)) {
    message(label)
  }
  state <- new.env(parent = emptyenv())
  state$bar <- utils::txtProgressBar(min = 0L, max = total, style = 3L)
  state$open <- TRUE
  step <- max(1L, as.integer(total) %/% 100L)
  list(
    tick = function(iteration) {
      if (state$open && (iteration %% step == 0L || iteration == total)) {
        utils::setTxtProgressBar(state$bar, iteration)
      }
      invisible(NULL)
    },
    close = function() {
      if (state$open) {
        state$open <- FALSE
        base::close(state$bar)
      }
      invisible(NULL)
    }
  )
}

# Progress-reporting wrapper for the raw-permutation engines, which re-index G
# directly instead of permuting Freedman-Lane residuals and are already O(n^2).
.permuted_statistics <- function(total, compute, label = NULL) {
  statistics <- numeric(total)
  if (!total) {
    return(statistics)
  }
  progress <- .permutation_progress(total, label)
  on.exit(progress$close(), add = TRUE)
  for (row in seq_len(total)) {
    statistics[[row]] <- compute(row)
    progress$tick(row)
  }
  statistics
}
