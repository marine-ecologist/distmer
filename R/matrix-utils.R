#' Validate and return a full distance matrix
#'
#' @param x A `dist` object or a numeric square matrix.
#' @param tol Numerical tolerance for symmetry and diagonal checks.
#' @return A numeric symmetric matrix.
#' @export
as_distance_matrix <- function(x, tol = sqrt(.Machine$double.eps)) {
  if (inherits(x, "dist")) {
    x <- as.matrix(x)
  }

  if (!is.matrix(x) || !is.numeric(x)) {
    stop("`x` must be a numeric square matrix or a `dist` object.", call. = FALSE)
  }

  if (nrow(x) != ncol(x)) {
    stop("Distance matrix must be square.", call. = FALSE)
  }

  if (nrow(x) < 2L) {
    stop("Distance matrix must contain at least two observations.", call. = FALSE)
  }

  if (any(!is.finite(x))) {
    stop("Distance matrix contains non-finite values.", call. = FALSE)
  }

  if (any(x < -tol)) {
    stop("Distance matrix contains negative values.", call. = FALSE)
  }

  if (max(abs(x - t(x))) > tol) {
    stop("Distance matrix is not symmetric within tolerance.", call. = FALSE)
  }

  if (max(abs(diag(x))) > tol) {
    stop("Distance matrix diagonal is not zero within tolerance.", call. = FALSE)
  }

  x[x < 0 & x >= -tol] <- 0
  x <- (x + t(x)) / 2
  diag(x) <- 0
  x
}

#' Diagnose the Euclidean spectrum of a dissimilarity
#'
#' Negative eigenvalues are retained and summarized rather than silently
#' discarded. Corrections, when added, will be separate sensitivity analyses.
#'
#' @param distance A `dist` object or full distance matrix.
#' @param tol Eigenvalue tolerance.
#' @return A list containing the centred matrix, complete eigenvalue spectrum,
#'   inertia summaries, and a Euclidean flag.
#' @export
distance_diagnostics <- function(distance, tol = 1e-10) {
  G <- gower_center(distance)
  eigenvalues <- eigen(G, symmetric = TRUE, only.values = TRUE)$values
  negative <- eigenvalues < -tol
  positive <- eigenvalues > tol
  positive_inertia <- sum(eigenvalues[positive])
  negative_inertia <- -sum(eigenvalues[negative])

  list(
    G = G,
    eigenvalues = eigenvalues,
    positive_inertia = positive_inertia,
    negative_inertia = negative_inertia,
    negative_count = sum(negative),
    negative_ratio = if (positive_inertia > 0) negative_inertia / positive_inertia else NA_real_,
    euclidean = !any(negative),
    correction = "none"
  )
}

#' Gower-centre a distance matrix
#'
#' Implements `G = -0.5 * J %*% D^2 %*% J`, where
#' `J = I - 11'/N`.
#'
#' @param distance A `dist` object or full distance matrix.
#' @return The centred inner-product matrix `G`.
#' @export
gower_center <- function(distance) {
  distance <- as_distance_matrix(distance)
  n <- nrow(distance)
  J <- diag(n) - matrix(1 / n, nrow = n, ncol = n)
  G <- -0.5 * J %*% (distance ^ 2) %*% J
  (G + t(G)) / 2
}

#' Orthogonal projection matrix
#'
#' Uses a rank-revealing QR decomposition and therefore does not require an
#' explicit matrix inverse.
#'
#' @param X Numeric design matrix.
#' @param tol QR tolerance.
#' @return A list containing the projection matrix and estimated rank.
#' @export
projection_matrix <- function(X, tol = 1e-10) {
  X <- as.matrix(X)

  if (!is.numeric(X)) {
    stop("`X` must be numeric.", call. = FALSE)
  }

  n <- nrow(X)

  if (ncol(X) == 0L) {
    return(list(H = matrix(0, n, n), rank = 0L))
  }

  qr_X <- qr(X, tol = tol, LAPACK = FALSE)
  rank <- qr_X$rank

  if (rank == 0L) {
    return(list(H = matrix(0, n, n), rank = 0L))
  }

  Q <- qr.Q(qr_X, complete = FALSE)[, seq_len(rank), drop = FALSE]
  H <- tcrossprod(Q)
  H <- (H + t(H)) / 2

  list(H = H, rank = rank)
}

# Project columns into a declared model stratum without allowing floating-point
# remnants of a structurally zero component to create a spurious rank. Columns
# are scaled after projection, so the result depends on their span rather than
# their measurement units.
.projection_in_subspace <- function(H, X, tol = 1e-10) {
  X <- as.matrix(X)
  n <- nrow(H)
  if (!is.numeric(X) || nrow(X) != n || ncol(H) != n) {
    stop("`H` and `X` must be conformable numeric matrices.", call. = FALSE)
  }
  if (!ncol(X)) {
    return(list(H = matrix(0, n, n), rank = 0L,
                projected = matrix(0, n, 0L), active = logical()))
  }

  projected <- H %*% X
  original_norm <- sqrt(colSums(X^2))
  projected_norm <- sqrt(colSums(projected^2))
  relative_floor <- max(100 * tol, sqrt(.Machine$double.eps))
  active <- original_norm > 0 &
    projected_norm > relative_floor * original_norm
  projected[, !active] <- 0
  if (!any(active)) {
    return(list(H = matrix(0, n, n), rank = 0L,
                projected = matrix(0, n, ncol(X)), active = active))
  }

  basis <- sweep(projected[, active, drop = FALSE], 2L,
                 projected_norm[active], "/")
  out <- projection_matrix(basis, tol = relative_floor)
  out$projected <- projected
  out$active <- active
  out
}

trace_product <- function(A, B) {
  sum(A * t(B))
}

trace_hgh <- function(H, G) {
  # For a symmetric idempotent H, tr(H G H) = tr(H G).
  trace_product(H, G)
}
