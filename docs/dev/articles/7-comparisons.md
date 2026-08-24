# 7. Validation and reproducibility

A distance-based mixed model should be reproducible at three levels: the
data and ecological distance, the experimental design and permutation
rule, and the software result. `distmer` exposes each of these rather
than treating a pseudo-F and p-value as a black box.

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

## 

| Study | Design | Agreement |
|:---|:---|:---|
| Lejeune et al. 2024 | one-way, two treatments, 36 response traits | 36 of 36 models exact |
| Mitchell et al. 2025 | crossed fixed factors, unequal replication | 27 of 27 terms exact |
| Turner et al. 2025 | three fixed factors, islands random, split-plot by depth | all df exact; 27 of 27 within-cell contrasts exact |
| Micaroni et al. 2021 | chambers in treatment, sponges in chambers, crossed with time | every row of both models exact |
| Pintar & Dorn 2025 | BACI with nested random sites and three covariates | all df and denominators exact; SS median 1.4% out, worst 15% |

### One-way: Lejeune et al. (2024)

Euclidean distance on a single trait, untransformed, 9999 permutations.

``` r

ponds <- read_delim("Macroinvertebrates-taxonomic-and-functional-data.csv",
                    delim = ";", show_col_types = FALSE) |>
  mutate(Treatment = factor(Treatment, levels = c("Newt", "Fish")))

D <- ponds |> pull(Predator) |> dist(method = "euclidean")

permanova(D ~ Treatment, data = ponds, permutations = 9999, seed = 1)
```

### Crossed fixed factors: Mitchell et al. (2025)

`log(x + 1)`, Euclidean, Type III partial sums of squares. Region has
four levels with unequal replication, so Type III and Type II differ
here — and Type III is what PRIMER reports.

``` r

density <- read_excel("MHW_Foraging_DATA.xlsx", sheet = "Australia Density") |>
  transmute(Region = factor(Region), Species = factor(Species),
            y = log(Density + 1))

D <- density |> pull(y) |> dist(method = "euclidean")

permanova(D ~ Region * Species, data = density,
          permutations = 9999, seed = 1, type = "III")
```

### Split-plot with a random unit: Turner et al. (2025)

Depth is measured within island, the other two factors between islands,
so the design carries three error strata. Square-root transform,
Euclidean distance.

``` r

benthos <- read_csv("Turner_et_al_2025_Benthos_2010_2014.csv",
                    show_col_types = FALSE) |>
  transmute(
    ISLAND     = factor(ISLAND),
    ECOREGION  = factor(ECOREGION),
    POP_STATUS = factor(POP_STATUS),
    DEPTH      = factor(DEPTH_BIN_NEW, levels = c("Shallow", "Mid", "Deep")),
    across(c(`CORAL_%`, `CCA_%`, `SC_%`, `HALIMEDA_%`, `MA_%`, `EMA_%`,
             `TURF_%`), sqrt)
  )

D <- benthos |> select(ends_with("%")) |> dist(method = "euclidean")

permanova(D ~ DEPTH * ECOREGION * POP_STATUS + (1 | ISLAND) + (1 | ISLAND:DEPTH),
          data = benthos, permutations = 999, seed = 1)
```

The published pairwise table for this study compares depths *within*
each ecoregion-by-population-status cell, which is the family described
in the contrasts vignette:

``` r

permanova_contrasts(fit, type = "pairwise", factor = "DEPTH",
                    within = c("ECOREGION", "POP_STATUS"), permutations = 9999)
```

### Nested repeated measures: Micaroni et al. (2021)

Chambers nested in treatment, sponges nested in chambers, both random,
crossed with time. Naming the nesting in the term itself —
`treatment:chamber` rather than `chamber` — is what lets the model see
that chamber sits inside treatment, and it is what separates a nested
random effect from a genuine interaction when the expected mean squares
are assembled.

``` r

sponges <- read_excel("Micaroni_et_al_sponge_hypoxia.xlsx", sheet = "2.13") |>
  mutate(treatment = factor(treatment, levels = c("100%a.s.", "5%a.s.")),
         chamber   = factor(chamber),
         sponge    = factor(sponge),
         time      = factor(time, levels = paste0("T", 0:12)))

# PRIMER reports Bray-Curtis as a percentage
D <- sponges |> transmute(y = sqrt(pumping_rate)) |> vegdist(method = "bray") * 100

permanova(D ~ treatment * time + (1 | treatment:chamber) +
            (1 | treatment:chamber:sponge) + (1 | treatment:chamber:time),
          data = sponges, permutations = 999, seed = 1)
```

Treatment is tested over chambers, chambers over sponges, and time and
the treatment-by-time interaction over the chamber-by-time stratum — all
derived from the design rather than declared.

### Covariates and a BACI contrast: Pintar & Dorn (2025)

Before/after crossed with impact/control, site-plots random within
region, and three hydrological covariates carried as slopes.

``` r

eels <- bind_rows(srs, tsl) |>
  mutate(across(where(is.numeric), \(x) replace_na(x, 0))) |>
  filter(PERIOD == 3, YEAR %in% c(1996:2001, 2015:2022)) |>
  mutate(REGION   = factor(REGION),
         SITEPLOT = factor(SITEPLOT),
         TS       = factor(if_else(YEAR <= 2001, "Before", "After"),
                           levels = c("Before", "After")))

D <- eels |> select(ADIXEN:PALPAL) |>
  vegdist(method = "jaccard", binary = TRUE) * 100

permanova_crossed_incomplete(
  compile_permanova_design(
    ~ DEPTH_AVE_30DAY + log.LDS + log.DSD + REGION * TS + (1 | REGION:SITEPLOT),
    eels, D),
  permutations = 999
)
```

### The expected-mean-square construction, without any data

PRIMER prints its own expected mean squares *and* the pseudo-F it builds
from them. That makes the denominator machinery testable on its own:
feed `.ems_f_construction()` PRIMER’s expected mean squares and check it
returns PRIMER’s pseudo-F and both degrees of freedom. Across the
workspaces above this covers **340 constructions** — fixed, crossed,
unbalanced mixed and nested repeated-measures designs, under Type I and
Type III — and all 340 reproduce PRIMER’s pseudo-F and numerator degrees
of freedom exactly.

Two of those constructions are pinned as fixtures in the package test
suite – the Turner split-plot and the Micaroni nested design, the latter
chosen because its sponge stratum and its residual have identical
expected mean squares, so only the nesting hierarchy can say which is
the denominator. Those run on every check because they carry PRIMER’s
numbers and need no data. The remaining constructions were swept against
the workspaces directly and are not part of the suite.

### What did not reproduce

Two designs could not be reproduced, in both cases because the published
data is not the data the workspace was built from:

- the male *Limnogonus* analysis uses 20 morphometric variables, two of
  which do not appear in the published measurements;
- two of the four *Micaroni* models rest on a sheet whose total sum of
  squares differs from the workspace by 0.03%.

The *Pintar & Dorn* reconstruction is partial for the same reason. The
workspace names its own 403 samples, so the design is recoverable
exactly – and the degrees of freedom and denominators above come out
exactly right – but a handful of rows still do not match the published
values, which leaves the total sum of squares 0.5% adrift and one
covariate 15% out.

Both are recorded rather than worked around. A benchmark that cannot be
reproduced from the archive is a fact about the archive.

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
