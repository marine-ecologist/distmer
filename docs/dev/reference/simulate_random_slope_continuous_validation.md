# Validate the continuous random-slope engine

[`permanova_random_slope_continuous()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_random_slope_continuous.md)
tests the equally weighted mean within-group slope of a numeric
covariate against the variation among group slopes. For a Euclidean
response that is exactly the one-sample test on the per-group ordinary
least squares slopes, so the reference here is `t.test(slopes)` squared.

## Usage

``` r
simulate_random_slope_continuous_validation(
  simulations = 200L,
  permutations = 99L,
  effect = 0.8,
  groups = 8L,
  occasions = 5L,
  slope_sd = 0.5,
  residual_sd = 0.3,
  alpha = 0.05,
  seed = 20260805L,
  keep_results = FALSE
)
```

## Arguments

- simulations:

  Iterations per scenario.

- permutations:

  Permutations per fit.

- effect:

  Location shift applied in the alternative scenario.

- groups:

  Number of random groups.

- occasions:

  Covariate values observed within every group.

- slope_sd:

  Standard deviation of the group slopes.

- residual_sd:

  Observation-level simulation standard deviation.

- alpha:

  Nominal level.

- seed:

  Base random seed.

- keep_results:

  Retain the per-iteration rows.

## Value

An object of class `permanova_random_slope_continuous_validation`.
