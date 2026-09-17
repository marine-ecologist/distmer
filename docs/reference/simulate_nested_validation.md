# Validate balanced nested mixed-model Type I error

Simulates treatment assigned to top-level regions with sites nested
within regions. The distance-based statistic is compared with ANOVA on
region means, and binomial uncertainty is reported for the null
rejection rate.

## Usage

``` r
simulate_nested_validation(
  simulations = 500L,
  regions_per_treatment = 4L,
  sites_per_region = 2L,
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

  Top-level units per treatment.

- sites_per_region:

  Nested sites per region.

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

A nested-validation report.
