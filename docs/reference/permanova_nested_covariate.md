# Covariate-adjusted balanced nested PERMANOVA

Tests categorical factors in a balanced `/`-generated nested hierarchy
while adjusting for one or more continuous covariates. Each factor is
tested for an effect *beyond* the covariate relationship: its numerator
is the sum of squares added over the covariate(s), and its denominator
is the same exact nested expected-mean-square stratum used by
[`permanova_nested()`](https://marine-ecologist.github.io/distmer/reference/permanova_nested.md).
The covariate is treated as an observation-level nuisance (charged to
the residual) and is itself tested over that residual.

## Usage

``` r
permanova_nested_covariate(
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

  Canonical Design Object whose fixed model contains at least one
  categorical factor and at least one continuous covariate, with one
  balanced `/`-generated random hierarchy.

- permutations:

  A non-negative integer or `"all"`.

- seed:

  Optional random seed.

- hypothesis:

  Sum-of-squares partition for the fixed factors: `"I"` (sequential),
  `"II"` (marginal), or `"III"` (partial). Types coincide for an
  orthogonal balanced fixed structure; Type III with interactions
  requires sum-to-zero contrasts, set automatically by the
  [`permanova()`](https://marine-ecologist.github.io/distmer/reference/permanova.md)
  front door.

- keep_permutations:

  Retain permutation statistics.

- max_exact:

  Maximum exact permutation space.

- tol:

  Numerical tolerance.

## Value

An object of class `permanova_nested`.

## Details

A single-slope covariate cannot be split cleanly across nested strata,
so a covariate with strong between-unit structure should be decomposed
by the user into explicit between- and within-unit predictors before
testing.
