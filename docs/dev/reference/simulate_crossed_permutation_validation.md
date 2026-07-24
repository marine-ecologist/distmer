# Validate crossed randomized-block permutation inference

Simulates a genuine treatment-by-block variance component under a null
fixed treatment effect. It compares the distance engine with a direct
response-level within-block cell permutation and the classical
interaction-F reference. A second scenario deliberately introduces
treatment-specific residual dispersion to expose exchangeability
sensitivity.

## Usage

``` r
simulate_crossed_permutation_validation(
  simulations = 500L,
  permutations = 99L,
  blocks = 8L,
  replicates = 3L,
  interaction_sd = 2,
  block_sd = 1,
  residual_sd = 1,
  dispersion_ratio = 2,
  alpha = 0.05,
  seed = 1L,
  keep_results = FALSE
)
```

## Arguments

- simulations:

  Number of simulations per scenario.

- permutations:

  Monte Carlo permutations per fitted model.

- blocks:

  Number of random blocks.

- replicates:

  Observations per treatment-by-block cell.

- interaction_sd, block_sd, residual_sd:

  Generating standard deviations.

- dispersion_ratio:

  Impact/control residual SD ratio in the sensitivity scenario.

- alpha:

  Rejection threshold.

- seed:

  Reproducibility seed.

- keep_results:

  Retain iteration-level results.

## Value

A crossed-permutation validation report.
