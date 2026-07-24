# Lower-level unbalanced nested PERMANOVA

Tests one categorical treatment assigned either to child units or to
observations in one unbalanced parent/child hierarchy. Child-level tests
use equally weighted child centroids, adjust for parent, and permute
intact children only within parent-by-exact-size strata.
Observation-level tests adjust for child and permute observation
residuals within child.

## Usage

``` r
permanova_unbalanced_nested_lower(
  design,
  level = c("auto", "child", "observation"),
  permutations = 999L,
  seed = NULL,
  keep_permutations = FALSE,
  max_exact = 100000L,
  tol = 1e-10
)
```

## Arguments

- design:

  Canonical unbalanced parent/child Design Object containing a distance.

- level:

  \`"auto"\`, \`"child"\`, or \`"observation"\`.

- permutations:

  A non-negative integer or \`"all"\`.

- seed:

  Optional random seed.

- keep_permutations:

  Retain permutation statistics.

- max_exact:

  Maximum exact permutation space.

- tol:

  Numerical tolerance.

## Value

An object of class \`permanova_unbalanced_nested_lower\`.
