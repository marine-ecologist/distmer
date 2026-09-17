# Estimate distance-based variance components

Uses the residual and each named random-effect ANOVA stratum to solve
the balanced trace-derived EMS system. Negative unrestricted estimates
are retained; truncation is an optional reported sensitivity analysis.

## Usage

``` r
estimate_variance_components(
  design,
  truncate = FALSE,
  assume_additive = FALSE,
  tol = 1e-10
)
```

## Arguments

- design:

  A balanced crossed canonical Design Object.

- truncate:

  Also report estimates truncated at zero.

- assume_additive:

  Explicitly assume omitted fixed-by-random interaction variances are
  zero. The default rejects such omissions.

- tol:

  Numerical tolerance.

## Value

A variance-component specification.
