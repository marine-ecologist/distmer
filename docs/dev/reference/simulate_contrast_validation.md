# Validate planned one-factor contrasts

Simulates unequal group counts and compares a planned A-minus-B distance
contrast with an independently constrained scalar model and the
classical one-df contrast F test under a null and location alternative.

## Usage

``` r
simulate_contrast_validation(
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

  B-versus-A location shift in the power scenario.

- residual_sd:

  Generating residual standard deviation.

- alpha:

  Rejection threshold.

- seed:

  Reproducibility seed.

- keep_results:

  Retain iteration-level results.

## Value

A Phase 13 validation report.
