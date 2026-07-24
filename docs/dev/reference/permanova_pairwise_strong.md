# Strong pairwise inference by intersection-specific closed permutation tests

Tests every pairwise hypothesis for a fixed one-factor design. For each
intersection null, observations are reassigned only within connected
components of the selected pairwise equalities; groups outside those
components remain fixed. A maximum pseudo-F tests that intersection, and
closure forms each adjusted p-value as the maximum intersection p-value
over all intersections containing that pair.

## Usage

``` r
permanova_pairwise_strong(
  design,
  permutations = 999L,
  seed = NULL,
  keep_permutations = FALSE,
  max_exact = 100000L,
  max_intersections = 1023L,
  tol = 1e-10
)
```

## Arguments

- design:

  Canonical fixed one-factor Design Object.

- permutations:

  A non-negative integer or \`"all"\`.

- seed:

  Optional random seed.

- keep_permutations:

  Retain intersection maximum statistics and label assignment orders.

- max_exact:

  Maximum unique label assignments for any intersection.

- max_intersections:

  Maximum intersection tests allowed by the computational guard.

- tol:

  Numerical tolerance.

## Value

An object of class \`permanova_pairwise_strong\`.

## Details

The finite-sample strong-FWER proof requires independent observations
and equality of the complete observation distributions within every
true-null connected component. It is not a proof for equality of
centroids alone when dispersions or other distributional features
differ. The method does not use or claim subset pivotality.
