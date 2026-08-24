# Validate the unbalanced split-plot mixed engine

Compares
[`permanova_unbalanced_splitplot()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_unbalanced_splitplot.md)
with an `lm` on equally weighted unit centroids for the whole-plot
factor, which is the estimand the engine declares under unequal
replication.

## Usage

``` r
simulate_unbalanced_splitplot_validation(
  simulations = 200L,
  permutations = 99L,
  effect = 1,
  unit_sd = 1,
  residual_sd = 0.5,
  alpha = 0.05,
  seed = 20260802L,
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

An object of class `permanova_unbalanced_splitplot_validation`.
