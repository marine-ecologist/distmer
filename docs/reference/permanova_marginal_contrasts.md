# Multifactor marginal contrasts

Tests contrasts among explicit marginal means from a full-rank
categorical fixed model. Each contrast has its own nuisance-retaining
reduced model and unrestricted Freedman–Lane residual permutation.

## Usage

``` r
permanova_marginal_contrasts(
  design,
  factor,
  contrasts,
  weights = c("equal", "proportional"),
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

  Canonical fixed-effect Design Object.

- factor:

  Target categorical predictor.

- contrasts:

  Numeric level contrast vector, row matrix, or named list.

- weights:

  `"equal"` or `"proportional"` nuisance-cell weights.

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

An object of class `permanova_marginal_contrasts`.
