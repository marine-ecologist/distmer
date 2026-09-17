# Benchmark balanced mixed inference

Benchmark balanced mixed inference

## Usage

``` r
benchmark_permanova(
  group_counts = c(8L, 16L),
  replicates = 3L,
  permutations = c(99L, 999L),
  seed = 1L
)
```

## Arguments

- group_counts:

  Whole-plot counts to benchmark (must be even).

- replicates:

  Observations per whole plot.

- permutations:

  Permutation counts to benchmark.

- seed:

  Reproducibility seed.

## Value

Timing data; results are machine-specific and make no performance claim.
