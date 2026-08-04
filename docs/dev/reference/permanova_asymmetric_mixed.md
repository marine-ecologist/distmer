# Asymmetric nested mixed PERMANOVA (unit + unit-by-subplot random strata)

Reproduces PRIMER's asymmetric nested mixed design with two random
terms: a unit nested in the crossed whole-plot factors, `(1 | unit)`,
and the unit by sub-plot-factor interaction, `(1 | unit:subplot)`. This
gives three error strata, all at the observation level (as in PRIMER,
which permutes residuals under a reduced model):

## Usage

``` r
permanova_asymmetric_mixed(
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

  A compiled Design Object with two random terms: one `(1 | unit)` main
  term and one `(1 | unit:subplot)` interaction term.

- permutations:

  A non-negative integer or `"all"`.

- seed:

  Optional random seed.

- hypothesis:

  Sum-of-squares partition for the fixed terms (`"I"`, `"II"`, `"III"`).

- keep_permutations:

  Retain permutation statistics on each test spec.

- max_exact:

  Maximum permutation space allowed for exact enumeration.

- tol:

  Numerical tolerance.

## Value

An object of class `permanova_nested` (shared table/summary shape).

## Details

- **IS(unit)** = between-unit variation after the whole-plot fixed terms
  – the denominator for the whole-plot fixed effects (those not
  involving the sub-plot factor);

- **subplot x IS(unit)** = the unit-by-sub-plot interaction – the
  denominator for the sub-plot fixed effects (those involving the
  sub-plot factor);

- **Residual** = within unit-by-sub-plot cells – the denominator for
  both random strata.

Degrees of freedom are derived from the actual (possibly incomplete)
design, so units missing a sub-plot level are handled. For the Caribbean
benthos example this is
`bray ~ DEPTH * ECOREGION * POP_STATUS + (1 | ISLAND) + (1 | ISLAND:DEPTH)`,
giving `IS(EC x PS)` and `DE x IS(EC x PS)`.
