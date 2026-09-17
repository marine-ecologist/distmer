# Validate unequal-replication random-intercept inference

Compares the Phase 11 distance engine with an independent scalar
cluster-centroid Freedman–Lane implementation under a null and location
alternative. Cluster sizes are 2, 3, and 5, with both treatments
represented within every size stratum.

## Usage

``` r
simulate_unequal_random_intercept_validation(
  simulations = 500L,
  permutations = 99L,
  effect = 1,
  cluster_sd = 1,
  residual_sd = 1,
  alpha = 0.05,
  seed = 20260719L,
  keep_results = FALSE
)
```

## Arguments

- simulations:

  Number of simulations per scenario.

- permutations:

  Monte Carlo permutations per fit.

- effect:

  Cluster-level treatment shift in the power scenario.

- cluster_sd, residual_sd:

  Generating standard deviations.

- alpha:

  Rejection threshold.

- seed:

  Reproducibility seed.

- keep_results:

  Retain iteration-level results.

## Value

A Phase 11 validation report.
