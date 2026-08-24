# Validate the covariate-adjusted balanced nested engine

[`permanova_nested_covariate()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_nested_covariate.md)
charges the covariate to the observation residual and tests each factor
for the sum of squares it adds *over* the covariate, against the exact
nested stratum of the factor-only model. This validator reconstructs
both pieces independently: the numerator as `SS(factor | covariate)`
from an observation-level `lm`, and the denominator as the corresponding
stratum mean square of `aov(y ~ tr + Error(reg/site))`.

## Usage

``` r
simulate_nested_covariate_validation(
  simulations = 200L,
  permutations = 99L,
  effect = 1,
  region_sd = 1,
  site_sd = 0.5,
  residual_sd = 0.4,
  slope = 0.7,
  alpha = 0.05,
  seed = 20260803L,
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

- region_sd, site_sd:

  Simulation standard deviations for the hierarchy.

- residual_sd:

  Observation-level simulation standard deviation.

- slope:

  Covariate slope used to generate the response.

- alpha:

  Nominal level.

- seed:

  Base random seed.

- keep_results:

  Retain the per-iteration rows.

## Value

An object of class `permanova_nested_covariate_validation`.
