# Pairwise one-factor distance comparisons

Constructs every pairwise difference among the levels of the sole fixed
factor and delegates to
[`permanova_contrasts()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_contrasts.md).

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

  Multiplicity method for
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) (default
  "holm").

- permutations:

  Number of permutations; defaults to `fit`'s value.

- seed:

  RNG seed; defaults to `fit`'s value.

- keep_permutations:

  Retain permutation statistics.

- max_exact:

  Maximum exact permutation space.

- tol:

  Numerical tolerance.

## Value

An object of class `permanova_contrasts` with a pairwise family.
