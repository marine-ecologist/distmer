# Contrasts for supported mixed PERMANOVA designs

Tests planned contrasts among a single categorical treatment's
experimental- unit centroids. The unit and exact interaction/residual
denominator are selected from the declared random structure. Supported
structures are one balanced or unequally replicated random intercept, a
top-level treatment in one balanced or unbalanced parent/child
hierarchy, a complete or incomplete treatment-by-block design with its
declared random interaction, and the balanced or validated unbalanced
treatment-by-region/site hybrid.

## Usage

``` r
permanova_mixed_contrasts(
  design,
  contrasts,
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

- contrasts:

  Numeric treatment-level contrast vector, row matrix, or named list.

- adjust:

  Multiplicity adjustment.

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

An object of class `permanova_mixed_contrasts`.
