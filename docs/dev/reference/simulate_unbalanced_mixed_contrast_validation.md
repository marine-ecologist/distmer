# Validate unequal-replication mixed planned contrasts

Simulates a three-treatment random-intercept design in which every exact
cluster-size stratum contains all treatments. The A-minus-B contrast is
compared with a scalar GLH on equally weighted cluster centroids and
with a direct reconstruction of the restricted Freedman–Lane
distribution.

## Usage

``` r
simulate_unbalanced_mixed_contrast_validation(
  simulations = 500L,
  permutations = 99L,
  effect = 0.8,
  alpha = 0.05,
  seed = 20260720L,
  keep_results = FALSE
)
```

## Arguments

- simulations:

  Number of simulations per scenario.

- permutations:

  Permutations per fit.

- effect:

  A-minus-B location shift in the power scenario.

- alpha:

  Rejection threshold.

- seed:

  Reproducibility seed.

- keep_results:

  Retain simulation-level results.

## Value

An unbalanced mixed-contrast validation report.
