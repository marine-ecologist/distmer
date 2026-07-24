# Validate balanced crossed-nested hybrid inference

Simulates a fixed treatment crossed with regions and sites nested within
regions. The engine is compared with ANOVA on treatment-by-region
profile means and an independent response-level implementation of every
restricted treatment-profile permutation.

## Usage

``` r
simulate_hybrid_validation(
  simulations = 500L,
  permutations = 99L,
  regions = 8L,
  sites = 3L,
  replicates = 2L,
  region_sd = 1,
  site_sd = 0.7,
  interaction_sd = 1.2,
  nested_interaction_sd = 0.6,
  residual_sd = 0.4,
  effect = 1,
  alpha = 0.05,
  seed = 20260718L,
  keep_results = FALSE
)
```

## Arguments

- simulations:

  Number of simulations per null and location scenario.

- permutations:

  Monte Carlo treatment-profile permutations per fit.

- regions, sites, replicates:

  Balanced hierarchy dimensions.

- region_sd, site_sd, interaction_sd, nested_interaction_sd,
  residual_sd:

  Generating standard deviations.

- effect:

  Fixed treatment shift in the power scenario.

- alpha:

  Rejection threshold.

- seed:

  Reproducibility seed.

- keep_results:

  Retain iteration-level results.

## Value

A hybrid-validation report.
