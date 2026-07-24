# Simulate a balanced random-intercept community design

Simulate a balanced random-intercept community design

## Usage

``` r
simulate_balanced_mixed(
  groups_per_treatment = 4L,
  replicates = 3L,
  response_variables = 1L,
  effect = 0,
  cluster_sd = 1,
  residual_sd = 1,
  dispersion_ratio = 1,
  seed = NULL
)
```

## Arguments

- groups_per_treatment:

  Whole plots per treatment level.

- replicates:

  Observations per whole plot.

- response_variables:

  Number of response dimensions.

- effect:

  Euclidean location-shift magnitude.

- cluster_sd:

  Random-intercept standard deviation.

- residual_sd:

  Residual standard deviation in the control group.

- dispersion_ratio:

  Residual SD multiplier in the impact group.

- seed:

  Optional random seed.

## Value

A list with model data, response matrix, and Euclidean distance.
