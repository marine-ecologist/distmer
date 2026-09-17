# Unbalanced multi-factor split-plot PERMANOVA

Tests several categorical factors in an unequally replicated
`(1 | unit)` random-intercept design. Whole-plot factors (constant
within each unit) are tested on **equally weighted unit centroids** over
the between-unit stratum – the unit-weighted estimand of
[`permanova_unequal_random_intercept()`](https://marine-ecologist.github.io/distmer/reference/permanova_unequal_random_intercept.md),
generalised to any number of fixed factors and sums-of-squares type –
with size-restricted whole-unit permutation. Within-unit (sub-plot)
factors, and any between-by-within interaction, are tested over the
within-unit residual with within-unit residual permutation. This is the
distance-based analogue of PRIMER's asymmetric nested design (for
example `IS(EC x PS)` as the whole-plot stratum and `DE x IS(EC x PS)`
as the sub-plot residual) when islands are unequally sampled.

## Usage

``` r
permanova_unbalanced_splitplot(
  design,
  permutations = 999L,
  seed = NULL,
  hypothesis = "I",
  keep_permutations = FALSE,
  max_exact = 100000L,
  tol = 1e-10
)
```

## Arguments

- design:

  A compiled Design Object with one `(1 | unit)` random term and one or
  more categorical fixed factors.

- permutations:

  A non-negative integer or `"all"`.

- seed:

  Optional random seed.

- hypothesis:

  Sum-of-squares partition for the fixed terms (`"I"`, `"II"`, `"III"`),
  applied within the whole-plot and within-unit blocks separately.

- keep_permutations:

  Retain permutation statistics on each test spec.

- max_exact:

  Maximum permutation space allowed for exact enumeration.

- tol:

  Numerical tolerance.

## Value

An object of class `permanova_nested` (shared table/summary shape).
