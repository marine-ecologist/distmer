# Balanced crossed-nested hybrid PERMANOVA

Tests one fixed treatment crossed with regions and with sites nested
within regions. Required covariance components are region, region:site,
treatment:region, and treatment:region:site. Complete
treatment-by-region profiles permute within regions and carry every
nested descendant intact.

## Usage

``` r
permanova_hybrid(
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

  Canonical Phase 9 hybrid Design Object.

- permutations:

  A non-negative integer or \`"all"\`.

- seed:

  Optional permutation seed.

- keep_permutations:

  Retain permutation statistics.

- truncate:

  Also report zero-truncated variance components.

- max_exact:

  Maximum exact permutation space.

- tol:

  Numerical tolerance.

## Value

An object of class \`permanova_hybrid\`.
