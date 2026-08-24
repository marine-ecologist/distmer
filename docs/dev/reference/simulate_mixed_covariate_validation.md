# Validate the balanced split-plot mixed engine

Compares
[`permanova_mixed_covariate()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_mixed_covariate.md)
with the classical split-plot `aov(y ~ tr * sub + Error(unit))`
decomposition: the whole-plot factor is tested in the unit stratum and
the sub-plot factor and interaction in the within-unit stratum.

## Usage

``` r
simulate_mixed_covariate_validation(
  simulations = 200L,
  permutations = 99L,
  effect = 1,
  unit_sd = 1,
  residual_sd = 0.5,
  alpha = 0.05,
  seed = 20260801L,
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

- unit_sd, residual_sd:

  Simulation standard deviations.

- alpha:

  Nominal level.

- seed:

  Base random seed.

- keep_results:

  Retain the per-iteration rows.

## Value

An object of class `permanova_mixed_covariate_validation`.
