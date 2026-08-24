# Multi-level unbalanced nested PERMANOVA

Tests one or more fixed terms in an unbalanced `/`-generated hierarchy
of any depth. The two-level engines
([`permanova_unbalanced_nested()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_unbalanced_nested.md)
and
[`permanova_unbalanced_nested_lower()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_unbalanced_nested_lower.md))
each test a single categorical treatment in one parent/child pair; this
engine covers deeper chains such as `(1 | region/site/quadrat)` and any
number of crossed fixed factors, interactions, and covariates.

## Usage

``` r
permanova_unbalanced_nested_multi(
  design,
  permutations = 999L,
  seed = NULL,
  hypothesis = "I",
  keep_permutations = FALSE,
  max_exact = 100000L,
  tol = 1e-10
)
```

## Arguments

- design:

  Compiled Design Object with one contiguous `/`-generated hierarchy of
  at least two levels.

- permutations:

  A non-negative integer or `"all"`.

- seed:

  Optional random seed.

- hypothesis:

  Sum-of-squares partition applied among the terms sharing a hierarchy
  level (`"I"`, `"II"`, or `"III"`).

- keep_permutations:

  Retain the permutation statistics on each test spec.

- max_exact:

  Maximum permutation space allowed for exact enumeration.

- tol:

  Numerical tolerance.

## Value

An object of class `permanova_nested`.

## Details

Each fixed term is assigned to the shallowest hierarchy level within
which it is constant, and tested on **equally weighted centroids at that
level** over the residual among those centroids after the parent level
and the other terms sharing the level. Aggregation cascades up the chain
– observations onto the deepest units, those onto their parents with
equal weight per child, and so on – so a level's centroid is not
reweighted by how many observations happen to sit beneath it. A term
that varies inside the deepest unit (typically an observation-level
covariate) is tested over the observation residual with residuals
permuted within the deepest units.

Permutation moves intact units within their own parent, and only across
units whose complete subtree shape matches, so an unequal hierarchy
never exchanges structurally different units. No unbalanced
expected-mean-square approximation is made.
