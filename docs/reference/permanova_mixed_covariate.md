# Split-plot / randomized-block / ANCOVA PERMANOVA for one random intercept

Tests any number of categorical factors, and optionally continuous
covariates, in a balanced `(1 | group)` random-intercept design. Each
fixed term is classified relative to the grouping unit and tested over
the matching stratum:

## Usage

``` r
permanova_mixed_covariate(
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

  Canonical Design Object with one `(1 | group)` random term and one or
  more fixed factors (a continuous covariate is optional).

- permutations:

  A non-negative integer or `"all"`.

- seed:

  Optional random seed.

- hypothesis:

  Sum-of-squares partition for the fixed terms: `"I"` (sequential),
  `"II"` (marginal, respecting containment), or `"III"` (partial). Types
  coincide for an orthogonal balanced fixed structure and differ only
  when fixed terms are non-orthogonal (interactions). Type III with
  interactions requires sum-to-zero contrasts, which the
  [`permanova()`](https://marine-ecologist.github.io/distmer/reference/permanova.md)
  front door sets automatically.

- keep_permutations:

  Retain permutation statistics on each test spec.

- max_exact:

  Maximum permutation space allowed for exact enumeration.

- tol:

  Numerical tolerance.

## Value

An object of class `permanova_nested` (shared table/summary shape).

## Details

- **whole-plot** – a factor constant within every unit is tested over
  the between-unit random-intercept stratum, exactly as in
  [`permanova_mixed()`](https://marine-ecologist.github.io/distmer/reference/permanova_mixed.md);

- **within-unit** – a sub-plot / randomized-block factor fully crossed
  within every unit, a covariate, or a covariate-by-factor slope
  interaction is tested over the within-unit residual;

- **between-by-within** – an interaction with any within-unit component
  is routed to the within-unit residual;

- **partial** – a factor neither constant nor fully crossed within the
  unit is genuinely confounded with the random unit and is refused, so
  no design receives a silently wrong denominator.

This is the distance-based split-plot analysis: whole-plot effects over
the unit stratum, sub-plot effects and covariates over the residual,
covering multi-factor split-plot, randomized-block, and
analysis-of-covariance layouts. A covariate that varies within units is
tested over the residual while its between-unit component is retained as
nuisance adjustment for whole-plot terms
([`decompose_covariate()`](https://marine-ecologist.github.io/distmer/reference/decompose_covariate.md)
separates between- and within-unit slopes when required). The
single-whole-plot-factor case without a covariate is handled by the
exact-EMS
[`permanova_mixed()`](https://marine-ecologist.github.io/distmer/reference/permanova_mixed.md);
the
[`permanova()`](https://marine-ecologist.github.io/distmer/reference/permanova.md)
front door dispatches automatically.
