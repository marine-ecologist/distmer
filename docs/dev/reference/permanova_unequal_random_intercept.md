# Unequal-replication random-intercept PERMANOVA

Tests one cluster-level categorical treatment using equally weighted
random- intercept unit centroids. Freedman–Lane residuals are permuted
only among whole clusters having the same observation count. The target
is therefore a cluster-weighted location effect; no observation-weighted
EMS or parametric F distribution is inferred.

## Usage

``` r
permanova_unequal_random_intercept(
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

An object of class \`permanova_unequal_random_intercept\`.
