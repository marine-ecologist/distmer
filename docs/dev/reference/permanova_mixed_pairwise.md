# Pairwise contrasts for supported mixed PERMANOVA designs

Pairwise contrasts for supported mixed PERMANOVA designs

## Usage

``` r
permanova_mixed_pairwise(
  design,
  adjust = c("holm", "bonferroni", "BH", "BY", "none"),
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

- adjust:

  Multiplicity adjustment.

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

An object of class \`permanova_mixed_contrasts\`.
