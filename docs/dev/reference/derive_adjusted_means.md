# Derive explicitly covariate-adjusted means

Constructs coefficient rows for a target factor while holding every
continuous nuisance variable at a user-supplied value. Categorical
nuisance factors are averaged either over their complete reference grid
or over observed joint frequencies. No covariate reference value is
inferred.

## Usage

``` r
derive_adjusted_means(
  design,
  factor,
  at = NULL,
  weights = c("equal", "proportional"),
  tol = 1e-10
)
```

## Arguments

- design:

  Canonical fixed-effect Design Object.

- factor:

  Target categorical predictor.

- at:

  Named numeric vector or list giving one finite scalar for every
  continuous nuisance variable.

- weights:

  `"equal"` or `"proportional"` categorical nuisance weights.

- tol:

  Numerical tolerance.

## Value

Adjusted-mean coefficient rows and reference-grid metadata.
