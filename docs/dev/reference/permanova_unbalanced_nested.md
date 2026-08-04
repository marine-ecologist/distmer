# Unbalanced nested top-unit PERMANOVA

Tests one treatment assigned to top-level units of one unbalanced nested
hierarchy. Equally weighted top-unit centroids are permuted only among
hierarchies with the same sorted vector of child observation counts.

## Usage

``` r
permanova_unbalanced_nested(
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

  Canonical Design Object with one random intercept.

- permutations:

  A non-negative integer or `"all"`.

- seed:

  Optional random seed.

- keep_permutations:

  Retain permutation statistics.

- max_exact:

  Maximum exact permutation space.

- tol:

  Numerical tolerance.

## Value

An object of class `permanova_unbalanced_nested`.
