# Pairwise multifactor marginal comparisons

Pairwise multifactor marginal comparisons

## Usage

``` r
permanova_marginal_pairwise(
  design,
  factor,
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

- weights:

  \`"equal"\` or \`"proportional"\` nuisance-cell weights.

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

An object of class \`permanova_marginal_contrasts\`.
