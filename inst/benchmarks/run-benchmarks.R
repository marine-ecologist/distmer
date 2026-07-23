library(distmer)

# Machine-specific results are printed, not committed as universal claims.
print(benchmark_permanova(
  group_counts = c(8L, 16L, 32L),
  replicates = 3L,
  permutations = c(99L, 999L),
  seed = 20260718L
))
