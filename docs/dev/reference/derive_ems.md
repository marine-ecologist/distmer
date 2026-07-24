# Derive EMS for a balanced random-intercept design

The supported Phase 2 envelope contains exactly one random intercept,
equal replication per grouping unit, and fixed predictors that are
constant within grouping units.

## Usage

``` r
derive_ems(design, hypothesis = "I", tol = 1e-10)
```

## Arguments

- design:

  A canonical \`permanova_design\` containing a distance matrix.

- hypothesis:

  Fixed-effect sum-of-squares partition: \`"I"\`, \`"II"\`, or
  \`"III"\`.

- tol:

  Numerical tolerance.

## Value

A named list of symbolic EMS specifications.
