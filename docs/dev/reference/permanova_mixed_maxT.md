# Global-null max-F inference for supported mixed designs

Uses one shared raw experimental-unit permutation distribution for a
family of treatment contrasts. Units and strata are selected from the
declared random structure. Adjusted p-values compare each observed
pseudo-F with the permutation maximum. This controls familywise error
under the complete treatment-null randomization hypothesis (weak FWER);
strong subset-null control is not claimed.

## Usage

``` r
permanova_mixed_maxT(
  design,
  contrasts,
  permutations = 999L,
  seed = NULL,
  keep_permutations = FALSE,
  max_exact = 100000L,
  tol = 1e-10
)
```

## Arguments

- design:

  Canonical mixed Design Object containing a distance.

- contrasts:

  Numeric treatment-level contrast vector, row matrix, or named list.

- permutations:

  A non-negative integer or \`"all"\`.

- seed:

  Optional random seed.

- keep_permutations:

  Retain the joint statistic matrix and maxima.

- max_exact:

  Maximum exact permutation space.

- tol:

  Numerical tolerance.

## Value

An object of class \`permanova_mixed_maxT\`.
