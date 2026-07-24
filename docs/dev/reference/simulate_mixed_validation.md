# Run mixed-model operating-characteristic simulations

Evaluates null rejection, location-shift power, dispersion sensitivity,
denominator agreement with scalar ANOVA, failures, runtime, and
numerical stability. Rejection-rate uncertainty is reported as a Monte
Carlo standard error; small runs are smoke tests, not inferential
validation.

## Usage

``` r
simulate_mixed_validation(
  simulations = 100L,
  permutations = 199L,
  alpha = 0.05,
  groups_per_treatment = 4L,
  replicates = 3L,
  effect = 1,
  dispersion_ratio = 2,
  seed = 1L,
  keep_results = FALSE
)
```

## Arguments

- simulations:

  Number of simulations per scenario.

- permutations:

  Permutations per fitted model.

- alpha:

  Rejection threshold.

- groups_per_treatment, replicates:

  Balanced design sizes.

- effect:

  Location shift for the power scenario.

- dispersion_ratio:

  Impact/control residual SD ratio for the dispersion scenario.

- seed:

  Simulation seed.

- keep_results:

  Retain simulation-level results.

## Value

A simulation-validation report.
