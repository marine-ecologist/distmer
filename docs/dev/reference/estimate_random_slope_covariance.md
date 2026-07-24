# Estimate a distance-based random intercept–slope covariance

Returns trace covariance components for the subject midpoint (random
intercept at centred condition coding) and the second-minus-first
condition difference (random slope). Equal-replication within-cell noise
is removed by method of moments. Negative unrestricted components are
retained; optional eigenvalue truncation is a separate sensitivity.

## Usage

``` r
estimate_random_slope_covariance(design, truncate = FALSE, tol = 1e-10)
```

## Arguments

- design:

  Canonical Phase 27 random-slope Design Object.

- truncate:

  Also return a positive-semidefinite eigenvalue-truncated sensitivity
  matrix.

- tol:

  Numerical tolerance.

## Value

Random-coefficient covariance estimates and diagnostics.
