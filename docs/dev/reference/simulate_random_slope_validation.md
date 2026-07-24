# Validate balanced two-condition random-slope inference

Compares the distance engine with paired-difference scalar F tests,
direct whole-subject sign-flip reconstruction, and scalar
method-of-moments random intercept–slope covariance estimates. A centred
asymmetric null is included as a labelled sensitivity outside the
central-symmetry guarantee.

## Usage

``` r
simulate_random_slope_validation(
  simulations = 500L,
  permutations = 99L,
  effect = 0.8,
  alpha = 0.05,
  seed = 20260725L,
  keep_results = FALSE
)
```

## Arguments

- simulations:

  Number of simulations per scenario.

- permutations:

  Sign flips per fit.

- effect:

  Fixed mean difference in the power scenario.

- alpha:

  Rejection threshold.

- seed:

  Reproducibility seed.

- keep_results:

  Retain simulation-level results.

## Value

A random-slope validation report.
