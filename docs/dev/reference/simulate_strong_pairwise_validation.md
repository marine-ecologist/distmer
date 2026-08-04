# Validate strong closed pairwise permutation inference

Audits complete- and partial-null familywise rejection, reconstructs
every sampled intersection maximum directly from the retained label
assignments, and verifies that the partial `A = B` intersection never
moves observations from the non-null `C` group.

## Usage

``` r
simulate_strong_pairwise_validation(
  simulations = 500L,
  permutations = 99L,
  effect = 1.2,
  alpha = 0.05,
  seed = 20260726L,
  keep_results = FALSE
)
```

## Arguments

- simulations:

  Number of simulations per scenario.

- permutations:

  Label assignments per intersection.

- effect:

  Location shift for the non-null group in the partial-null scenario.

- alpha:

  Familywise rejection threshold.

- seed:

  Reproducibility seed.

- keep_results:

  Retain simulation-level results.

## Value

A strong-pairwise validation report.
