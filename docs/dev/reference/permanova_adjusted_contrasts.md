# Explicitly covariate-adjusted distance contrasts

Tests contrasts among adjusted means with contrast-specific
Freedman–Lane residual permutations and the common full-model residual
denominator.

## Usage

``` r
permanova_adjusted_contrasts(
  design,
  factor,
  contrasts,
  at = NULL,
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

- at:

  Named numeric vector or list giving one finite scalar for every
  continuous nuisance variable.

- weights:

  \`"equal"\` or \`"proportional"\` categorical nuisance weights.

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

An object of class \`permanova_adjusted_contrasts\`.
