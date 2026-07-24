# Pairwise one-factor distance comparisons

Constructs every pairwise difference among the levels of the sole fixed
factor and delegates to \`permanova_contrasts()\`.

## Usage

``` r
permanova_pairwise(
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

  Canonical fixed one-factor Design Object.

- adjust:

  Multiplicity method: \`"holm"\`, \`"bonferroni"\`, \`"BH"\`, \`"BY"\`,
  or \`"none"\`.

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

An object of class \`permanova_contrasts\` with a pairwise family.
