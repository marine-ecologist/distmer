# Validate explicitly covariate-adjusted contrasts

Simulates an unbalanced fixed design with a categorical nuisance factor
and a continuous covariate. An adjusted A-minus-B contrast at an
explicit covariate value is compared with a scalar Gaussian GLH and with
a direct reconstruction of every sampled Freedman–Lane statistic.

## Usage

``` r
simulate_adjusted_contrast_validation(
  simulations = 500L,
  permutations = 99L,
  effect = 0.8,
  alpha = 0.05,
  at = 0.25,
  seed = 20260721L,
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

- at:

  Continuous-covariate reference value.

- seed:

  Reproducibility seed.

- keep_results:

  Retain simulation-level results.

## Value

A covariate-adjusted contrast validation report.
