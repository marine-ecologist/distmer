# Audit unequal-random-intercept calibration across independent batches

Uses exact unique treatment assignments within cluster-size strata and a
separate Monte Carlo approximation. Batch seeds are separated by one
million to prevent the overlapping per-iteration streams that make
adjacent-seed summaries appear artificially stable. Gaussian,
heavy-tail, unequal- allocation, and deliberate dispersion-violation
scenarios are included.

## Usage

``` r
audit_unequal_random_intercept_calibration(
  batches = 20L,
  simulations_per_batch = 500L,
  permutations = 99L,
  engine_checks = 10L,
  cluster_sd = 1,
  residual_sd = 1,
  alpha = 0.05,
  seed = 20260719L,
  keep_results = FALSE
)
```

## Arguments

- batches:

  Number of independent seed batches.

- simulations_per_batch:

  Null simulations in each batch and scenario.

- permutations:

  Monte Carlo assignment draws per simulation.

- engine_checks:

  Number of first-batch Gaussian datasets also passed through the
  complete distance engine with exhaustive permutations.

- cluster_sd, residual_sd:

  Generating standard deviations.

- alpha:

  Rejection threshold.

- seed:

  Base reproducibility seed. Batch seeds are spaced internally.

- keep_results:

  Retain simulation-level p-values.

## Value

A Phase 14 exact-randomization calibration audit.
