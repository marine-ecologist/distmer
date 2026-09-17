# Validate unbalanced crossed-nested hybrid inference

Simulates one treatment crossed with an unequal region/site hierarchy
with one missing treatment-by-region profile and unequal replicated site
cells. The distance engine is compared with the scalar equal-site-weight
profile model and direct reconstruction of every sampled profile
permutation.

## Usage

``` r
simulate_unbalanced_hybrid_validation(
  simulations = 500L,
  permutations = 99L,
  effect = 0.8,
  alpha = 0.05,
  seed = 20260724L,
  keep_results = FALSE
)
```

## Arguments

- simulations:

  Number of simulations per scenario.

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

An unbalanced hybrid validation report.
