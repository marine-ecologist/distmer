# Unbalanced crossed-nested hybrid PERMANOVA

Tests one treatment crossed with a region/site hierarchy when regions
have unequal site counts, treatment-by-region profiles are missing, or
replicated treatment-by-site cells have unequal sizes. Each observed
profile must still contain every site in its region. A profile is the
equal-weight mean of its site-cell centroids; observed profiles are
adjusted for region and permuted only within region and identical
site-specific cell-count topology.

## Usage

``` r
permanova_unbalanced_hybrid(
  design,
  permutations = 999L,
  seed = NULL,
  keep_permutations = FALSE,
  max_exact = 100000L,
  tol = 1e-10
)
```

## Arguments

- design:

  Canonical unbalanced hybrid Design Object with a distance.

- permutations:

  A non-negative integer or `"all"`.

- seed:

  Optional permutation seed.

- keep_permutations:

  Retain permutation statistics.

- max_exact:

  Maximum exact permutation space.

- tol:

  Numerical tolerance.

## Value

An object of class `permanova_unbalanced_hybrid`.
