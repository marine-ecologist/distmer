# Validate mixed-design global-null max-F inference

Simulates an unequally replicated three-treatment random-intercept
design whose exact-size strata each contain every treatment. The
complete shared statistic matrix is reconstructed from scalar unit
centroids, and family rejection is reported for unadjusted, Holm, and
weak-FWER max-F procedures.

## Usage

``` r
simulate_mixed_maxT_validation(
  simulations = 500L,
  permutations = 99L,
  effect = 0.8,
  alpha = 0.05,
  seed = 20260722L,
  keep_results = FALSE
)
```

## Arguments

- simulations:

  Number of simulations per scenario.

- permutations:

  Permutations per fit.

- effect:

  Group-C location shift in the power scenario.

- alpha:

  Family rejection threshold.

- seed:

  Reproducibility seed.

- keep_results:

  Retain simulation-level results.

## Value

A mixed global-null max-F validation report.
