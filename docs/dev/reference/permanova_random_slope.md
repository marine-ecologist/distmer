# Balanced two-condition random-slope PERMANOVA

Tests the population mean within-subject condition difference for a
correlated \`(1 + condition \| subject)\` model. Complete difference
vectors are sign-flipped by subject. Exactness is conditional on central
symmetry of subject difference vectors under the fixed-effect null; this
is stated in every result and is not presented as unrestricted label
exchangeability.

## Usage

``` r
permanova_random_slope(
  design,
  permutations = 999L,
  seed = NULL,
  keep_permutations = FALSE,
  truncate = FALSE,
  max_exact = 100000L,
  tol = 1e-10
)
```

## Arguments

- design:

  Canonical Phase 27 random-slope Design Object.

- permutations:

  A non-negative integer or \`"all"\`.

- seed:

  Optional sign-flip seed.

- keep_permutations:

  Retain sign-flip statistics.

- truncate:

  Also return a positive-semidefinite covariance sensitivity.

- max_exact:

  Maximum exact sign-flip space.

- tol:

  Numerical tolerance.

## Value

An object of class \`permanova_random_slope\`.
