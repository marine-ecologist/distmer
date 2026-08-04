# Fixed-effect sequential PERMANOVA

Initial clean-room multi-factor engine. It implements Type-I/sequential
sums of squares and tests each term over the final residual mean square.
Permutations are unrestricted raw-data permutations represented by the
same row-and-column permutation of the Gower-centred matrix.

## Usage

``` r
permanova_fixed(
  distance,
  formula,
  data,
  permutations = 999L,
  seed = NULL,
  keep_permutations = FALSE,
  tol = 1e-10
)
```

## Arguments

- distance:

  A `dist` object or full distance matrix.

- formula:

  One-sided or two-sided model formula; the response is ignored.

- data:

  Model data in the same row order as the distance matrix.

- permutations:

  Number of random raw-data permutations.

- seed:

  Optional random seed.

- keep_permutations:

  Retain the matrix of permutation pseudo-F values.

- tol:

  QR tolerance.

## Value

An object of class `permanova_fixed`.

## Details

This is not yet a replacement for the legacy program's EMS denominator
logic for mixed, random, or nested designs.
