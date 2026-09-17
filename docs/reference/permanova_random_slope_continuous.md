# Continuous-covariate random-slope PERMANOVA

Tests the population slope of a continuous covariate that also carries a
random slope across groups, i.e. `distance ~ x + (1 + x | group)` with
`x` numeric – the covariate-by-random model of Anderson, Gorley & Clarke
(2008, Fig. 1.48b).

## Usage

``` r
permanova_random_slope_continuous(
  design,
  permutations = 999L,
  seed = NULL,
  keep_permutations = FALSE,
  max_exact = 100000L,
  tol = 1e-10
)
```

## Arguments

- design:

  A compiled Design Object with exactly one `(1 + x | group)` random
  term whose slope variable is numeric.

- permutations:

  A non-negative integer or `"all"` (whole-group sign flips).

- seed:

  Optional sign-flip seed.

- keep_permutations:

  Retain the sign-flip statistics.

- max_exact:

  Maximum exact sign-flip space.

- tol:

  Numerical tolerance.

## Value

An object of class `permanova_random_slope_continuous`.

## Details

The fixed slope has no exact single-source pseudo-F denominator in the
full observation space, because its expected mean square mixes the
random-slope and residual components with different coefficients. The
test therefore aggregates, exactly as the two-level categorical engine
reduces to per-subject difference vectors: it forms one within-group
slope vector per group (in the distance embedding), then tests whether
the equally weighted mean slope vector differs from zero against the
between-group variation of those slope vectors. Inference is by
whole-group sign flips of the slope vectors, valid under central
symmetry of the group slopes under the null. A parametric `F(1, ng - 1)`
reference is added only for rank-one Euclidean slope geometry with equal
within-group covariate information (`SSx`), where the statistic is the
classical one-sample squared t statistic under normal, identically
distributed group slopes. Euclidean geometry alone is not enough to
justify that reference for multivariate or heteroscedastic slopes.

The within-group slope is invariant to centring `x` within a group, so
the raw or a within-centred covariate
([`decompose_covariate()`](https://marine-ecologist.github.io/distmer/reference/decompose_covariate.md))
give the same test; the covariate need only vary inside each group.
