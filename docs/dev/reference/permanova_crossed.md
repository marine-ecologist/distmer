# Balanced crossed distance-based mixed model

Derives numerical EMS coefficients from covariance and hypothesis
projections, chooses exact mean-square denominators where possible, and
otherwise constructs an EMS-matched quasi-F denominator.

## Usage

``` r
permanova_crossed(
  design,
  p_method = c("none", "satterthwaite", "permutation"),
  truncate = FALSE,
  assume_additive = FALSE,
  tol = 1e-10,
  permutations = 999L,
  seed = NULL,
  keep_permutations = FALSE,
  max_exact = 100000L
)
```

## Arguments

- design:

  Canonical complete balanced crossed-factor design.

- p_method:

  `"none"` or opt-in `"satterthwaite"` F approximation.

- truncate:

  Also report zero-truncated variance-component estimates.

- assume_additive:

  Explicitly assume omitted fixed-by-random interaction variances are
  zero. The default rejects such omissions.

- tol:

  Numerical tolerance.

- permutations:

  A non-negative integer or `"all"` when `p_method = "permutation"`.

- seed:

  Optional permutation seed.

- keep_permutations:

  Retain permutation statistics.

- max_exact:

  Maximum exact permutation space.

## Value

An object of class `permanova_crossed`.
