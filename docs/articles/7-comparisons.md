# 7. Validation and reproducibility

``` r
knitr::opts_chunk$set(
  collapse = TRUE, comment = "#>", fig.width = 6.5, fig.height = 4.5,
  fig.align = "center",
  # the published-data chunks are slow (Turner alone is several minutes), so
  # they are cached beside the article; delete this folder to force a re-run,
  # e.g. after a change to distmer that would alter their output
  cache.path = "7-comparisons_cache/"
)
library(distmer)
library(tidyverse)
library(vegan)
library(ggplot2) 
library(kableExtra)
library(readxl)
theme_set(theme_minimal(base_size = 12))
```

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

## Published studies:

### One-way: Lejeune et al. (2024)

Euclidean distance on a single trait, untransformed, 9999 permutations.

``` r
ponds <- read_delim("~/distmer/code/distmer_tests/Lejeune_et_al/Macroinvertebrates-taxonomic-and-functional-data.csv",
                    delim = ";", show_col_types = FALSE) |>
  mutate(Treatment = factor(Treatment, levels = c("Newt", "Fish")))

D <- ponds |> pull(Predator) |> dist(method = "euclidean")

permanova(D ~ Treatment, data = ponds, permutations = 9999, seed = 1)
#> Distance-based PERMANOVA
#> Engine: fixed effects  (Type I SS)
#> 
#> Factors
#>  Name      Type  Levels
#>  Treatment Fixed 2     
#> 
#> PERMANOVA table of results
#>  Source    df SS    MS    Pseudo-F P(perm) P(MC)  Unique perms
#>  Treatment 1  0.071 0.071 3.5200   0.0688  0.0678 9999        
#>  Residual  30 0.604 0.020                                     
#>  Total     31 0.675
```

### Crossed fixed factors: Mitchell et al. (2025)

`log(x + 1)`, Euclidean, Type III partial sums of squares. Region has
four levels with unequal replication, so Type III and Type II differ
here — and Type III is what PRIMER reports.

``` r
density <- read_excel("~/distmer/code/distmer_tests/Mitchell_et_al/MHW_Foraging_DATA.xlsx", sheet = "Australia Density") |>
  transmute(Region = factor(Region), Species = factor(Species),
            y = log(Density + 1))
#> New names:
#> • `` -> `...4`

D <- density |> pull(y) |> dist(method = "euclidean")

permanova(D ~ Region * Species, data = density,
          permutations = 9999, seed = 1, type = "III")
#> Distance-based PERMANOVA
#> Engine: fixed effects  (Type III SS)
#> 
#> Factors
#>  Name           Type  Levels
#>  Region         Fixed  4    
#>  Species        Fixed  2    
#>  Region:Species Fixed NA    
#> 
#> PERMANOVA table of results
#>  Source         df  SS      MS    Pseudo-F P(perm) P(MC)  Unique perms
#>  Region         3   7.022   2.341 4.6929   0.0031  0.0033 9999        
#>  Species        1   2.558   2.558 5.1281   0.0260  0.0250 9999        
#>  Region:Species 3   2.045   0.682 1.3665   0.2595  0.2597 9999        
#>  Residual       204 101.744 0.499                                     
#>  Total          211 115.763
```

### Split-plot with a random unit: Turner et al. (2025)

Depth is measured within island, the other two factors between islands,
so the design carries three error strata. Square-root transform,
Euclidean distance.

``` r
benthos <- read_csv("~/distmer/code/distmer_tests/Turner_et_al/Turner_et_al_2025_Benthos_2010_2014.csv",
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
#> Distance-based PERMANOVA
#> Engine: asymmetric nested mixed  (Type I SS)
#> 
#> Factors
#>  Name                       Type                Levels
#>  DEPTH                      Fixed                3    
#>  ECOREGION                  Fixed                6    
#>  POP_STATUS                 Fixed                2    
#>  DEPTH:ECOREGION            Fixed               NA    
#>  DEPTH:POP_STATUS           Fixed               NA    
#>  ECOREGION:POP_STATUS       Fixed               NA    
#>  DEPTH:ECOREGION:POP_STATUS Fixed               NA    
#>  ISLAND                     Random, main        33    
#>  ISLAND:DEPTH               Random, interaction 97    
#> 
#> PERMANOVA table of results
#>  Source                     df   SS        MS      Pseudo-F P(perm)
#>  DEPTH                      2    435.187   224.757 6.7624   0.0010 
#>  ECOREGION                  5    226.786   45.357  14.3466  0.0010 
#>  POP_STATUS                 1    5.117     5.117   1.6186   0.1910 
#>  DEPTH:ECOREGION            10   653.683   67.852  2.9651   0.0010 
#>  DEPTH:POP_STATUS           2    501.897   254.868 9.8345   0.0010 
#>  ECOREGION:POP_STATUS       2    21.213    10.606  3.3549   0.0290 
#>  DEPTH:ECOREGION:POP_STATUS 4    204.508   51.127  3.4354   0.0020 
#>  IS(ISLAND)                 24   3072.784  128.033 13.1264         
#>  DEPTH x IS(ISLAND)         46   787.202   17.113  1.7545          
#>  Residual                   2142 20892.752 9.754                   
#>  Total                      2238 40867.347                         
#>  Unique perms
#>  999         
#>  999         
#>  999         
#>  999         
#>  999         
#>  999         
#>  999         
#>  999         
#>  999         
#>              
#> 
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
sponges <- read_excel("~/distmer/code/distmer_tests/Micaroni_et_al/Micaroni_et_al_sponge_hypoxia.xlsx", sheet = "2.13") |>
  mutate(treatment = factor(treatment, levels = c("100%a.s.", "5%a.s.")),
         chamber   = factor(chamber),
         sponge    = factor(sponge),
         time      = factor(time, levels = paste0("T", 0:12)))

# PRIMER reports Bray-Curtis as a percentage
D <- sponges |> transmute(y = sqrt(pumping_rate)) |> vegdist(method = "bray") * 100

permanova(D ~ treatment * time + (1 | treatment:chamber) +
            (1 | treatment:chamber:sponge) + (1 | treatment:chamber:time),
          data = sponges, permutations = 999, seed = 1)
#> Distance-based PERMANOVA
#> Engine: crossed mixed (incomplete / unequally replicated)
#> 
#> Factors
#>  Name                     Type                Levels
#>  treatment                Fixed                2    
#>  time                     Fixed               13    
#>  treatment:time           Fixed               NA    
#>  treatment:chamber        Random, interaction  6    
#>  treatment:chamber:sponge Random, interaction 18    
#>  treatment:chamber:time   Random, interaction 78    
#> 
#> PERMANOVA table of results
#>  Source                   df  SS         MS        Pseudo-F P(perm)
#>  treatment                1   42555.875  42555.875 50.8720  0.0010 
#>  time                     12  78718.177  6559.848  13.4035  0.0010 
#>  treatment:time           12  39319.903  3276.659  6.6951   0.0010 
#>  treatment:chamber        4   3346.114   836.529   0.1760   0.9990 
#>  treatment:chamber:sponge 12  57025.605  4752.134  8.0350   0.0010 
#>  treatment:chamber:time   48  23491.833  489.413   0.8275   0.9010 
#>  Residual                 144 85165.846  591.429                   
#>  Total                    233 329623.355                           
#>  Unique perms
#>              
#>              
#>              
#>              
#>              
#>              
#>              
#> 
```

Treatment is tested over chambers, chambers over sponges, and time and
the treatment-by-time interaction over the chamber-by-time stratum — all
derived from the design rather than declared.

## Run the deterministic checks

[`run_permanova_validation()`](https://marine-ecologist.github.io/distmer/reference/run_permanova_validation.md)
performs the package’s fast analytic battery. It checks reference
statistics, exact permutation counts, projector algebra, retained
non-Euclidean geometry, mixed-design denominators, and failure
boundaries.

``` r
# validation_report <- run_permanova_validation()
# validation_report[, c("check", "passed", "metric")]
# stopifnot(all(validation_report$passed))
```
