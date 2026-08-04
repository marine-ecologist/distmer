# The shared Freedman-Lane kernel replaces an explicit `A %*% G %*% t(A)` in
# every residual-permutation engine, so it is pinned directly against that
# definition rather than only through the engines that call it.

.fl_fixture <- function(n = 24L, groups = 4L, seed = 42L) {
  set.seed(seed)
  Y <- matrix(rnorm(n * 7L), n, 7L)
  G <- tcrossprod(Y) - 0.3 * diag(n)          # retains negative eigenvalues
  G <- (G + t(G)) / 2
  group <- factor(rep_len(seq_len(groups), n))
  full <- projection_matrix(stats::model.matrix(~ group))$H
  reduced <- projection_matrix(matrix(1, n, 1L))$H
  list(
    n = n, G = G, group = group, H_full = full, H_reduced = reduced,
    H_numerator = full - reduced, H_denominator = diag(n) - full,
    numerator_df = groups - 1L, denominator_df = n - groups
  )
}

# The definition the kernel is a fast path for.
.fl_reference <- function(fixture, index, H) {
  A <- fixture$H_reduced +
    (diag(fixture$n) - fixture$H_reduced)[index, , drop = FALSE]
  sum(H * t(A %*% fixture$G %*% t(A)))
}

test_that("kernel reproduces A G A' for low- and high-rank projectors", {
  fx <- .fl_fixture()
  kernel <- .freedman_lane_kernel(
    fx$G, fx$H_reduced, list(fx$H_numerator, fx$H_denominator)
  )
  set.seed(1)
  for (i in seq_len(25L)) {
    index <- sample.int(fx$n)
    expect_equal(
      kernel$evaluate(index),
      c(.fl_reference(fx, index, fx$H_numerator),
        .fl_reference(fx, index, fx$H_denominator)),
      tolerance = 1e-9
    )
  }
})

test_that("kernel reproduces A G A' under restricted within-block permutation", {
  fx <- .fl_fixture()
  blocks <- split(seq_len(fx$n), fx$group)
  kernel <- .freedman_lane_kernel(
    fx$G, fx$H_reduced, list(fx$H_numerator, fx$H_denominator)
  )
  set.seed(2)
  for (i in seq_len(25L)) {
    index <- seq_len(fx$n)
    for (block in blocks) index[block] <- sample(block)
    expect_equal(
      kernel$evaluate(index),
      c(.fl_reference(fx, index, fx$H_numerator),
        .fl_reference(fx, index, fx$H_denominator)),
      tolerance = 1e-9
    )
  }
})

test_that("low-rank representation picks the narrower of H and its complement", {
  fx <- .fl_fixture()
  numerator <- .low_rank_representation(fx$H_numerator)
  denominator <- .low_rank_representation(fx$H_denominator)

  expect_equal(numerator$scale, 0)
  expect_equal(ncol(numerator$W), fx$numerator_df)
  # A residual projector is carried through its complement, not its own rank.
  expect_equal(denominator$scale, 1)
  expect_equal(ncol(denominator$W), fx$n - fx$denominator_df)

  expect_equal(
    numerator$W %*% (numerator$signs * t(numerator$W)),
    fx$H_numerator, tolerance = 1e-9
  )
  expect_equal(
    diag(fx$n) + denominator$W %*% (denominator$signs * t(denominator$W)),
    fx$H_denominator, tolerance = 1e-9
  )
})

test_that("kernel falls back to the dense path and stays exact", {
  fx <- .fl_fixture()
  # Neither H nor I - H is narrow, so no low-rank form exists.
  basis <- eigen(fx$H_numerator, symmetric = TRUE)$vectors[
    , seq_len(fx$n %/% 2L + 2L), drop = FALSE
  ]
  wide <- 2 * tcrossprod(basis)
  asymmetric <- matrix(stats::rnorm(fx$n * fx$n), fx$n, fx$n)

  expect_null(.low_rank_representation(wide))
  expect_null(.low_rank_representation(asymmetric))

  # Mixing dense and low-rank projectors in one kernel must not cross wires.
  kernel <- .freedman_lane_kernel(
    fx$G, fx$H_reduced, list(fx$H_numerator, wide, asymmetric)
  )
  set.seed(3)
  for (i in seq_len(15L)) {
    index <- sample.int(fx$n)
    expect_equal(
      kernel$evaluate(index),
      c(.fl_reference(fx, index, fx$H_numerator),
        .fl_reference(fx, index, wide),
        .fl_reference(fx, index, asymmetric)),
      tolerance = 1e-9
    )
  }
})

test_that("kernel handles a rank-zero projector", {
  fx <- .fl_fixture()
  kernel <- .freedman_lane_kernel(
    fx$G, fx$H_reduced, list(matrix(0, fx$n, fx$n))
  )
  expect_equal(kernel$evaluate(sample.int(fx$n)), 0)
})

test_that("statistics driver matches the explicit permutation loop", {
  fx <- .fl_fixture()
  set.seed(4)
  orders <- t(vapply(seq_len(50L), function(i) sample.int(fx$n), integer(fx$n)))

  explicit <- vapply(seq_len(nrow(orders)), function(row) {
    index <- orders[row, ]
    (.fl_reference(fx, index, fx$H_numerator) / fx$numerator_df) /
      (.fl_reference(fx, index, fx$H_denominator) / fx$denominator_df)
  }, numeric(1))

  expect_equal(
    .freedman_lane_statistics(
      fx$G, fx$H_reduced, fx$H_numerator, fx$numerator_df,
      fx$H_denominator, fx$denominator_df, orders
    ),
    explicit, tolerance = 1e-9
  )
})

test_that("statistics driver applies the block-expansion index", {
  fx <- .fl_fixture()
  blocks <- unname(split(seq_len(fx$n), fx$group))
  set.seed(5)
  orders <- t(vapply(seq_len(20L), function(i) sample.int(length(blocks)),
                     integer(length(blocks))))
  expand <- function(order) .block_permutation_index(blocks, order)

  explicit <- vapply(seq_len(nrow(orders)), function(row) {
    index <- expand(orders[row, ])
    (.fl_reference(fx, index, fx$H_numerator) / fx$numerator_df) /
      (.fl_reference(fx, index, fx$H_denominator) / fx$denominator_df)
  }, numeric(1))

  expect_equal(
    .freedman_lane_statistics(
      fx$G, fx$H_reduced, fx$H_numerator, fx$numerator_df,
      fx$H_denominator, fx$denominator_df, orders, index = expand
    ),
    explicit, tolerance = 1e-9
  )
})

test_that("an empty permutation plan returns no statistics", {
  fx <- .fl_fixture()
  expect_length(
    .freedman_lane_statistics(
      fx$G, fx$H_reduced, fx$H_numerator, fx$numerator_df,
      fx$H_denominator, fx$denominator_df,
      matrix(integer(), nrow = 0L, ncol = fx$n)
    ),
    0L
  )
})

test_that("progress reporting is opt-in and threshold-gated", {
  emitted <- function(total) {
    # The label is emitted at construction, so capture that too.
    output <- utils::capture.output({
      messages <- utils::capture.output({
        bar <- .permutation_progress(total, "label")
        for (i in seq_len(total)) bar$tick(i)
        bar$close()
      }, type = "message")
    }, type = "output")
    length(output) + length(messages) > 0L
  }

  withr_option <- function(value, code) {
    old <- getOption("distmer.progress")
    on.exit(options(distmer.progress = old), add = TRUE)
    options(distmer.progress = value)
    code
  }

  expect_false(withr_option(FALSE, emitted(500L)))
  expect_true(withr_option(TRUE, emitted(500L)))
  # Short runs stay silent even when reporting is switched on.
  expect_false(withr_option(TRUE, emitted(50L)))
  # Closing twice must not error or emit a second time.
  expect_silent({
    bar <- .permutation_progress(10L, NULL)
    bar$close()
    bar$close()
  })
})
