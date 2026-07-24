# Validate unbalanced nested and randomized-block inference

Simulates the two Phase 12 envelopes under a null and location
alternative, comparing every distance statistic and p-value with an
independent scalar top-unit or cell-profile reconstruction.

## Usage

``` r
simulate_unbalanced_mixed_validation(
  simulations = 500L,
  permutations = 99L,
  effect = 1,
  parent_sd = 1,
  child_sd = 0.7,
  block_sd = 1,
  interaction_sd = 1,
  residual_sd = 0.5,
  alpha = 0.05,
  seed = 20260719L,
  keep_results = FALSE
)
```

## Arguments

- simulations:

  Number of simulations per structure and scenario.

- permutations:

  Monte Carlo permutations per fit.

- effect:

  Treatment shift in the power scenario.

- parent_sd, child_sd, block_sd, interaction_sd, residual_sd:

  Generating standard deviations.

- alpha:

  Rejection threshold.

- seed:

  Reproducibility seed.

- keep_results:

  Retain iteration-level results.

## Value

A Phase 12 validation report.
