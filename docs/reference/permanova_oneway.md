# One-way PERMANOVA

Reconstructs the one-way algorithm in Anderson (2001). The `pairwise`
engine directly evaluates the published squared-distance equations and
therefore requires equal replication. The `projection` engine evaluates
the equivalent Gower-centred projection-matrix formulation and also
works with unequal group sizes.

## Usage

``` r
permanova_oneway(
  distance,
  group,
  permutations = 999L,
  seed = NULL,
  engine = c("pairwise", "projection"),
  keep_permutations = FALSE,
  tol = 1e-10
)
```

## Arguments

- distance:

  A `dist` object or full distance matrix.

- group:

  Grouping factor.

- permutations:

  Number of random permutations. Use `0` for no test.

- seed:

  Optional integer seed.

- engine:

  Either `"pairwise"` or `"projection"`.

- keep_permutations:

  Retain the permutation pseudo-F values.

- tol:

  Numerical tolerance for the numerically-zero residual check.

## Value

An object of class `permanova_oneway`.
