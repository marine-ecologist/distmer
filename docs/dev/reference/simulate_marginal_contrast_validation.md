# Validate multifactor marginal contrasts

Simulates an unequal-cell two-factor design and compares an equal-weight
marginal A contrast with a direct scalar general linear hypothesis and
response-level Freedman–Lane reconstruction.

## Usage

``` r
simulate_marginal_contrast_validation(
  simulations = 500L,
  permutations = 99L,
  effect = 0.8,
  residual_sd = 1,
  alpha = 0.05,
  seed = 20260719L,
  keep_results = FALSE
)
```

## Arguments

- simulations:

  Number of simulations per scenario.

- permutations:

  Monte Carlo permutations per fit.

- effect:

  Equal-grid marginal A shift in the power scenario.

- residual_sd:

  Generating residual standard deviation.

- alpha:

  Rejection threshold.

- seed:

  Reproducibility seed.

- keep_results:

  Retain iteration-level results.

## Value

A Phase 15 validation report.
