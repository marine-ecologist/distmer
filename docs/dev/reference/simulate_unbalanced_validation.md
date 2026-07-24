# Validate unbalanced fixed-effect inference

Simulates a fixed unequal-cell additive design and compares the Type II
distance engine with nested scalar linear models and an independent
response-level Freedman–Lane implementation under the null and a
location alternative.

## Usage

``` r
simulate_unbalanced_validation(
  simulations = 500L,
  permutations = 99L,
  effect = 0.8,
  residual_sd = 1,
  alpha = 0.05,
  seed = 20260718L,
  keep_results = FALSE
)
```

## Arguments

- simulations:

  Number of simulations per scenario.

- permutations:

  Monte Carlo residual permutations per fit.

- effect:

  Fixed A-factor effect in the power scenario.

- residual_sd:

  Generating residual standard deviation.

- alpha:

  Rejection threshold.

- seed:

  Reproducibility seed.

- keep_results:

  Retain iteration-level results.

## Value

An unbalanced-validation report.
