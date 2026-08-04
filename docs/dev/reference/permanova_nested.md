# Balanced nested distance-based mixed model

Assigns each fixed term to one orthogonal hierarchy stratum, uses that
stratum's exact residual mean square, and permutes units within their
parent while moving all descendants intact.

## Usage

``` r
permanova_nested(
  design,
  permutations = 999L,
  seed = NULL,
  keep_permutations = FALSE,
  truncate = FALSE,
  max_exact = 100000L,
  tol = 1e-10
)
```

## Arguments

- design:

  Canonical Design Object compiled with `distance`.

- permutations:

  A non-negative integer or `"all"`.

- seed:

  Optional random seed.

- keep_permutations:

  Retain permutation statistics.

- truncate:

  Also report zero-truncated nested variance components.

- max_exact:

  Maximum exact permutation space.

- tol:

  Numerical tolerance.

## Value

An object of class `permanova_nested`.
