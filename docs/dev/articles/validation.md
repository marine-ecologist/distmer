# 7. Validation and reproducibility

A distance-based mixed model should be reproducible at three levels: the
data and ecological distance, the experimental design and permutation
rule, and the software result. `distmer` exposes each of these rather
than treating a pseudo-F and p-value as a black box.

This guide shows reference comparisons, deterministic package checks,
simulation tools, and a reporting record that can accompany a paper or
data release.

## Reference agreement on vegan data

For an independent one-factor design, the fixed-effect statistic has the
same dissimilarity-space definition as
[`vegan::adonis2()`](https://vegandevs.github.io/vegan/reference/adonis.html).
The following comparison uses the `dune` communities and Bray–Curtis
dissimilarity.

``` r
data(dune, dune.env, package = "vegan")
dune_bray <- vegdist(dune, method = "bray")

distmer_fit <- permanova(
  dune_bray ~ Management, dune.env,
  permutations = 199, mc_permutations = 0, seed = 101
)
vegan_fit <- adonis2(
  dune_bray ~ Management, dune.env,
  permutations = 199
)

data.frame(
  implementation = c("distmer", "vegan::adonis2"),
  pseudo_f = c(distmer_fit$table$`Pseudo-F`[1], vegan_fit$F[1])
)  |> tab()
```

| implementation | pseudo_f |
|:---------------|---------:|
| distmer        |    2.767 |
| vegan::adonis2 |    2.767 |

The statistic should agree to numerical tolerance. Permutation p-values
need not be identical unless the two functions receive the identical
permutation matrix; recording only a seed does not guarantee that
different packages draw the same sequence of arrangements.

The dispersion companion can also be checked against
[`vegan::betadisper()`](https://vegandevs.github.io/vegan/reference/betadisper.html).

``` r
selected <- dune.env$Management %in% c("NM", "SF")
two_group_distance <- vegdist(dune[selected, ], method = "bray")
two_group_env <- droplevels(dune.env[selected, , drop = FALSE])
dispersion_design <- compile_permanova_design(
  ~ Management, two_group_env, two_group_distance
)
distmer_dispersion <- permanova_dispersion(
  dispersion_design, term = "Management", type = "median",
  permutations = 199, seed = 101
)
vegan_dispersion <- betadisper(
  two_group_distance, two_group_env$Management, type = "median"
)

c(
  maximum_distance_difference = max(abs(
    distmer_dispersion$distances - vegan_dispersion$distances
  )),
  pseudo_f_difference = distmer_dispersion$test$statistic -
    unname(anova(vegan_dispersion)["Groups", "F value"])
) |> tab()
```

|                             |   x |
|:----------------------------|----:|
| maximum_distance_difference |   0 |
| pseudo_f_difference         |   0 |

These are useful anchors, but validation is broader than agreement with
one package. Mixed designs are also checked against their univariate
Euclidean counterparts, exact enumeration, independent matrix
calculations, and operating-characteristic simulations.

## validation with PERMANOVA

## Run the deterministic checks

[`run_permanova_validation()`](https://marine-ecologist.github.io/distmer/dev/reference/run_permanova_validation.md)
performs the package’s fast analytic battery. It checks reference
statistics, exact permutation counts, projector algebra, retained
non-Euclidean geometry, mixed-design denominators, and failure
boundaries.

``` r
# validation_report <- run_permanova_validation()
# validation_report[, c("check", "passed", "metric")]
# stopifnot(all(validation_report$passed))
```

A green report shows that the installed code reproduces its pinned
numeric checks on the current platform; it does not certify that a new
field design satisfies the package assumptions.
