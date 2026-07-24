# Validate and return a full distance matrix

Validate and return a full distance matrix

## Usage

``` r
as_distance_matrix(x, tol = sqrt(.Machine$double.eps))
```

## Arguments

- x:

  A \`dist\` object or a numeric square matrix.

- tol:

  Numerical tolerance for symmetry and diagonal checks.

## Value

A numeric symmetric matrix.
