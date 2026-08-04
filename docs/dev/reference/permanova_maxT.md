# Global-null max-F contrast inference

Uses one shared raw-observation permutation distribution for every
planned contrast in a fixed one-factor design. Adjusted p-values compare
each observed pseudo-F with the permutation maximum. This controls
familywise error under the complete equal-centroid null (weak FWER); it
does not claim strong subset-null control.

## Usage

``` r
permanova_maxT(
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

  Canonical fixed one-factor Design Object.

- contrasts:

  Numeric contrast vector, row matrix, or named list.

- permutations:

  A non-negative integer or `"all"`.

- seed:

  Optional random seed.

- keep_permutations:

  Retain the joint statistic matrix and maxima.

- max_exact:

  Maximum exact permutation space.

- tol:

  Numerical tolerance.

## Value

An object of class `permanova_maxT`.
