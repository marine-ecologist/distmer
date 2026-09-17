# Validate level-specific balanced nested inference

Simulates null fixed effects assigned to regions, sites within regions,
and observations within sites. Each distance-based pseudo-F and its
reference F tail probability are compared with the corresponding
classical nested ANOVA calculation.

## Usage

``` r
simulate_nested_level_validation(
  simulations = 500L,
  regions_per_treatment = 4L,
  replicates = 3L,
  region_sd = 1.5,
  site_sd = 1,
  residual_sd = 0.5,
  alpha = 0.05,
  seed = 1L,
  keep_results = FALSE
)
```

## Arguments

- simulations:

  Number of simulation replicates.

- regions_per_treatment:

  Regions per top-level treatment.

- replicates:

  Observations per site.

- region_sd, site_sd, residual_sd:

  Generating standard deviations.

- alpha:

  Rejection threshold.

- seed:

  Reproducibility seed.

- keep_results:

  Retain iteration-level results.

## Value

A level-specific nested-validation report.
