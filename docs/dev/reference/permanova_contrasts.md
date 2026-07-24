# Planned one-factor distance contrasts

Tests one or more contrasts among the centroids of a single categorical
fixed factor. Each row supplies coefficients applied to unweighted group
centroids and must sum to zero. A contrast-specific reduced model
retains all orthogonal group information; Freedman–Lane residual
permutations are unrestricted. Holm adjustment is the default and
remains valid for non-orthogonal contrast families.

## Usage

``` r
permanova_contrasts(
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

  Canonical fixed one-factor Design Object.

- contrasts:

  Numeric contrast vector, matrix with contrasts in rows, or named list
  of numeric vectors.

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

An object of class \`permanova_contrasts\`.
