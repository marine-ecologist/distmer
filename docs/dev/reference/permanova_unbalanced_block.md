# Unbalanced randomized-block cell-profile PERMANOVA

Tests one treatment in a randomized block design with incomplete
treatment incidence and/or unequal replicated cell sizes.
Treatment:block cell centroids receive equal weight and move only within
block and exact cell-size strata. The residual treatment:block mean
square is descriptive; inference is empirical.

## Usage

``` r
permanova_unbalanced_block(
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

An object of class \`permanova_unbalanced_block\`.
