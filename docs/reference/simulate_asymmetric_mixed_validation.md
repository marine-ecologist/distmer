# Validate the asymmetric nested mixed engine

Compares
[`permanova_asymmetric_mixed()`](https://marine-ecologist.github.io/distmer/reference/permanova_asymmetric_mixed.md)
with `aov(y ~ EC * DE + Error(ISLAND / DE))`: the whole-plot factor is
tested in the unit stratum, and the sub-plot factor and its interaction
in the unit-by-sub-plot stratum.

## Usage

``` r
simulate_asymmetric_mixed_validation(
  simulations = 200L,
  permutations = 99L,
  effect = 1,
  island_sd = 1,
  island_depth_sd = 0.6,
  residual_sd = 0.4,
  alpha = 0.05,
  seed = 20260804L,
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

- island_sd, island_depth_sd:

  Simulation standard deviations for the two random strata.

- residual_sd:

  Observation-level simulation standard deviation.

- alpha:

  Nominal level.

- seed:

  Base random seed.

- keep_results:

  Retain the per-iteration rows.

## Value

An object of class `permanova_asymmetric_mixed_validation`.
