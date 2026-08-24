.oneway_pairwise_stat <- function(distance, group) {
  distance <- as_distance_matrix(distance)
  group <- droplevels(factor(group))

  n_total <- nrow(distance)
  n_groups <- nlevels(group)
  group_sizes <- table(group)

  if (length(group) != n_total) {
    stop("`group` length must equal the number of observations.", call. = FALSE)
  }

  if (n_groups < 2L) {
    stop("`group` must contain at least two levels.", call. = FALSE)
  }

  if (length(unique(as.integer(group_sizes))) != 1L) {
    stop(
      "The direct Anderson (2001) pairwise-distance equation requires a balanced one-way design.",
      call. = FALSE
    )
  }

  n_rep <- as.integer(group_sizes[[1]])
  d2 <- distance ^ 2
  ss_total <- sum(d2[upper.tri(d2)]) / n_total

  ss_within <- sum(vapply(
    levels(group),
    function(level) {
      idx <- which(group == level)
      block <- d2[idx, idx, drop = FALSE]
      sum(block[upper.tri(block)]) / n_rep
    },
    numeric(1)
  ))

  ss_among <- ss_total - ss_within
  df_among <- n_groups - 1L
  df_within <- n_total - n_groups
  ms_among <- ss_among / df_among
  ms_within <- ss_within / df_within
  pseudo_f <- ms_among / ms_within

  c(
    ss_total = ss_total,
    ss_among = ss_among,
    ss_within = ss_within,
    df_among = df_among,
    df_within = df_within,
    ms_among = ms_among,
    ms_within = ms_within,
    pseudo_f = pseudo_f
  )
}

.oneway_projection_stat <- function(distance, group) {
  distance <- as_distance_matrix(distance)
  group <- droplevels(factor(group))
  n <- nrow(distance)

  if (length(group) != n) {
    stop("`group` length must equal the number of observations.", call. = FALSE)
  }

  G <- gower_center(distance)
  X0 <- matrix(1, nrow = n, ncol = 1)
  X1 <- model.matrix(~ group)

  P0 <- projection_matrix(X0)
  P1 <- projection_matrix(X1)

  H_term <- P1$H - P0$H
  H_residual <- diag(n) - P1$H

  ss_among <- trace_hgh(H_term, G)
  ss_within <- trace_hgh(H_residual, G)
  ss_total <- sum(diag(G))
  df_among <- P1$rank - P0$rank
  df_within <- n - P1$rank
  ms_among <- ss_among / df_among
  ms_within <- ss_within / df_within
  pseudo_f <- ms_among / ms_within

  c(
    ss_total = ss_total,
    ss_among = ss_among,
    ss_within = ss_within,
    df_among = df_among,
    df_within = df_within,
    ms_among = ms_among,
    ms_within = ms_within,
    pseudo_f = pseudo_f
  )
}

#' One-way PERMANOVA
#'
#' Reconstructs the one-way algorithm in Anderson (2001). The `pairwise`
#' engine directly evaluates the published squared-distance equations and
#' therefore requires equal replication. The `projection` engine evaluates
#' the equivalent Gower-centred projection-matrix formulation and also works
#' with unequal group sizes.
#'
#' @param distance A `dist` object or full distance matrix.
#' @param group Grouping factor.
#' @param permutations Number of random permutations. Use `0` for no test.
#' @param seed Optional integer seed.
#' @param engine Either `"pairwise"` or `"projection"`.
#' @param keep_permutations Retain the permutation pseudo-F values.
#' @param tol Numerical tolerance for the numerically-zero residual check.
#' @return An object of class `permanova_oneway`.
#' @export
permanova_oneway <- function(
    distance,
    group,
    permutations = 999L,
    seed = NULL,
    engine = c("pairwise", "projection"),
    keep_permutations = FALSE,
    tol = 1e-10) {
  engine <- match.arg(engine)
  distance <- as_distance_matrix(distance)
  group <- droplevels(factor(group))

  if (length(group) != nrow(distance)) {
    stop("`group` length must equal the number of observations.", call. = FALSE)
  }

  if (length(permutations) != 1L || permutations < 0 || permutations != as.integer(permutations)) {
    stop("`permutations` must be one non-negative integer.", call. = FALSE)
  }

  stat_fun <- switch(
    engine,
    pairwise = .oneway_pairwise_stat,
    projection = .oneway_projection_stat
  )

  observed <- stat_fun(distance, group)
  # Match the mixed and unbalanced engines, which refuse rather than divide by a
  # numerically zero residual. Only the observed statistic is guarded: a single
  # degenerate permutation must not abort an otherwise valid run.
  if (.mean_square_is_zero(
    observed[["ms_within"]], tol = tol,
    reference = observed[["ss_total"]] / max(1L, nrow(distance) - 1L)
  )) {
    stop("The within-group residual mean square is numerically zero.",
         call. = FALSE)
  }
  permutation_f <- numeric(as.integer(permutations))

  if (permutations > 0L) {
    if (!is.null(seed)) {
      set.seed(seed)
    }

    permutation_f <- .permuted_statistics(
      as.integer(permutations),
      function(b) stat_fun(distance, sample(group, replace = FALSE))[["pseudo_f"]],
      label = "Permuting group labels"
    )

    p_value <- (1 + sum(permutation_f >= observed[["pseudo_f"]])) /
      (as.integer(permutations) + 1)
  } else {
    p_value <- NA_real_
  }

  table <- data.frame(
    term = c("group", "Residual", "Total"),
    df = c(observed[["df_among"]], observed[["df_within"]], nrow(distance) - 1L),
    ss = c(observed[["ss_among"]], observed[["ss_within"]], observed[["ss_total"]]),
    ms = c(observed[["ms_among"]], observed[["ms_within"]], NA_real_),
    pseudo_f = c(observed[["pseudo_f"]], NA_real_, NA_real_),
    p_perm = c(p_value, NA_real_, NA_real_),
    row.names = NULL
  )

  out <- list(
    call = match.call(),
    engine = engine,
    table = table,
    permutations = as.integer(permutations),
    seed = seed,
    unique_permutation_statistics = if (permutations > 0L) length(unique(permutation_f)) else 0L,
    permutation_f = if (keep_permutations) permutation_f else NULL
  )
  class(out) <- "permanova_oneway"
  out
}

#' @export
print.permanova_oneway <- function(x, ...) {
  cat("One-way PERMANOVA (", x$engine, " engine)\n", sep = "")
  print(x$table, row.names = FALSE)
  invisible(x)
}
