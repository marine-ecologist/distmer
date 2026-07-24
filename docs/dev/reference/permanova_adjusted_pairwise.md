# Pairwise explicitly covariate-adjusted distance contrasts

Pairwise explicitly covariate-adjusted distance contrasts

## Usage

``` r
permanova_adjusted_pairwise(
  design,
  factor,
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
