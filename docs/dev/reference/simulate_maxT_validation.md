# Validate global-null max-F multiplicity control

Simulates an unequal-size, three-group one-factor design under the
complete equal-centroid null and under a location shift. It reports
family rejection for unadjusted pairwise tests, Holm adjustment, and the
shared-permutation max-F adjustment.

## Usage

``` r
simulate_maxT_validation(
  simulations = 500L,
  permutations = 99L,
  effect = 0.8,
  alpha = 0.05,
  seed = 20260719L,
  keep_results = FALSE
)
```

## Arguments

- simulations:

  Number of simulations per scenario.

- permutations:

  Permutations per fit.

- effect:

  Location shift for group C in the power scenario.

- alpha:

  Rejection threshold.

- seed:

  Reproducibility seed.

- keep_results:

  Retain simulation-level results.

## Value

A global-null max-F validation report.
