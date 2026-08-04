# Validate design-aware multivariate dispersion inference

Evaluates equal-dispersion null calibration, invariance to a pure
location shift, and power against a dispersion ratio across one-way,
random-intercept, nested, and canonical randomized-block designs. Every
observed statistic is compared with an independently fitted scalar
classical reference.

## Usage

``` r
simulate_dispersion_validation(
  simulations = 200L,
  permutations = 99L,
  alpha = 0.05,
  seed = 20260718L,
  type = c("median", "centroid"),
  keep_results = FALSE
)
```

## Arguments

- simulations:

  Number of simulations per structure and scenario.

- permutations:

  Monte Carlo residual permutations per fit.

- alpha:

  Rejection threshold.

- seed:

  Reproducibility seed.

- type:

  Group spatial center: `"median"` or `"centroid"`.

- keep_results:

  Retain iteration-level results.

## Value

A dispersion-validation report.
