# Validate fixed-by-random interaction Type I error

Simulates a null treatment effect with genuine treatment-by-block random
variation. It compares the declared-interaction engine, the identical
classical interaction-denominator test, and an explicitly opted-in
additive model that omits the interaction variance.

## Usage

``` r
simulate_crossed_interaction_validation(
  simulations = 500L,
  blocks = 8L,
  replicates = 3L,
  interaction_sd = 2,
  block_sd = 1,
  residual_sd = 1,
  alpha = 0.05,
  seed = 1L,
  keep_results = FALSE
)
```

## Arguments

- simulations:

  Number of simulation replicates.

- blocks:

  Number of random blocks.

- replicates:

  Observations per treatment-by-block cell.

- interaction_sd, block_sd, residual_sd:

  Generating standard deviations.

- alpha:

  Rejection threshold.

- seed:

  Reproducibility seed.

- keep_results:

  Retain iteration-level results.

## Value

An interaction-validation report with binomial confidence intervals.
