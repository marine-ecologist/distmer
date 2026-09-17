# Validate mixed-design planned contrasts

Simulates a three-treatment randomized-block experiment with a genuine
treatment-by-block component. The planned A-minus-B contrast is checked
against a scalar general-linear-hypothesis test on treatment-by-block
cell means under both the null and a location alternative.

## Usage

``` r
simulate_mixed_contrast_validation(
  simulations = 500L,
  permutations = 99L,
  blocks = 8L,
  replicates = 2L,
  effect = 0.8,
  alpha = 0.05,
  seed = 20260719L,
  keep_results = FALSE
)
```

## Arguments

- simulations:

  Number of simulations per scenario.

- permutations:

  Permutations per fit.

- blocks:

  Number of random blocks.

- replicates:

  Observations per treatment-by-block cell.

- effect:

  A-minus-B location shift in the power scenario.

- alpha:

  Rejection threshold.

- seed:

  Reproducibility seed.

- keep_results:

  Retain simulation-level results.

## Value

A mixed-contrast validation report.
