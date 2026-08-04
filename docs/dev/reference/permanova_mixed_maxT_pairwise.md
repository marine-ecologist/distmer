# Pairwise global-null max-F inference for supported mixed designs

Pairwise global-null max-F inference for supported mixed designs

## Usage

``` r
permanova_mixed_maxT_pairwise(
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

  Canonical mixed Design Object containing a distance.

- permutations:

  A non-negative integer or `"all"`.

- seed:

  Optional random seed.

- keep_permutations:

  Retain the joint statistic matrix and maxima.

- max_exact:

  Maximum exact permutation space.

- tol:

  Numerical tolerance.

## Value

An object of class `permanova_mixed_maxT`.
