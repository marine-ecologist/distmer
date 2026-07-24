# Validate lower-level unbalanced nested inference

Simulates child-level and observation-level treatments in one unbalanced
parent/child hierarchy. Each distance result is compared with a scalar
unit- level general linear model and with direct reconstruction of the
restricted Freedman–Lane permutation distribution.

## Usage

``` r
simulate_unbalanced_nested_lower_validation(
  simulations = 500L,
  permutations = 99L,
  effect = 0.8,
  alpha = 0.05,
  seed = 20260723L,
  keep_results = FALSE
)
```

## Arguments

- simulations:

  Number of simulations per structure and scenario.

- permutations:

  Permutations per fit.

- effect:

  Treatment location shift in the power scenario.

- alpha:

  Rejection threshold.

- seed:

  Reproducibility seed.

- keep_results:

  Retain simulation-level results.

## Value

A lower-level unbalanced nested validation report.
