# Unbalanced fixed-effect PERMANOVA

Tests explicit estimable Type I, II, or III hypotheses over the final
full-model residual mean square. Type III requires a full-rank model and
explicit sum-to-zero contrasts when interactions are present. Type II
respects marginality: a term is adjusted for every term that does not
contain it and the higher-order interactions containing it are excluded,
so main effects remain estimable in the presence of an interaction
(matching the SAS/`car` "Type II" hypothesis).

## Usage

``` r
permanova_unbalanced(
  design,
  hypothesis = c("I", "II", "III"),
  terms = NULL,
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

- hypothesis:

  One of `"I"`, `"II"`, or `"III"`.

- terms:

  Optional fixed terms to test.

- permutations:

  A non-negative integer or `"all"`.

- seed:

  Optional permutation seed.

- keep_permutations:

  Retain permutation statistics.

- max_exact:

  Maximum exact permutation space.

- tol:

  Numerical tolerance.

## Value

An object of class `permanova_unbalanced`.
