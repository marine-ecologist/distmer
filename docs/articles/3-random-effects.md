# 3. Random effects and ecological sampling designs

Random effects describe ecological replication: `(1 | site)` says that
observations from the same site share a site-level community component.
Nested, blocked, and random-slope terms extend that idea to more
structured sampling.

- `(1 | site)` — repeated observations within sites; treatments assigned
  among sites are tested with sites as the independent units.
- `(1 | region/site)` — sites nested within regions; each predictor is
  tested at the hierarchy level where it varies.
- `(1 | block) + (1 | treatment:block)` — a treatment randomized within
  blocks; the declared treatment-by-block variation supplies the
  treatment denominator.
- `(1 | region/site) + (1 | treatment:region) + (1 | treatment:region:site)`
  — treatment crossed with a nested region/site hierarchy.
- `(1 + condition | subject)` or `(1 + x | group)` — correlated random
  intercepts and slopes for complete repeated units.

`distmer` deliberately borrows `lmer` formula notation because it is a
familiar and concise way to state this hierarchy. `distmer`applies the
declared fixed and random subspaces to a Gower-centred dissimilarity
matrix and obtains inference from design-appropriate denominators and
restricted permutations.

### Simple random effects:

Herbivory on coral reefs has a strong impact on macrobenthic community
succession. This can be demonstrated by using tiles placed in herbivore
exclusion cages, typically “open”, “partial” and “caged” treatments (see
[Evensen et al
2021](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.3536)
for details).

Below is changes in [community
succession](https://github.com/nickevensen/Palau-benthic-community-succession/)
through time using the open surfaces of tiles (“crowns”) across
different caging treatments:

``` r
library(tidyverse)
library(distmer)
library(patchwork)
library(ggrepel)


coralsuccession <- read.csv("~/distmer/code/Macro_communities.csv") |>
  filter(Nutrients == "Normal") |>
  filter(Microhabitat == "Crown")

env <- coralsuccession[,1:5]
coral <- coralsuccession[,6:17]
coral_bray <- vegdist(decostand(coral, method = "hellinger"), method = "bray")

coral_pcoa <- cmdscale(coral_bray, k = 2, eig = TRUE, add = TRUE)

coral_pcoa_scores <- as.data.frame(coral_pcoa$points) |>
  setNames(c("PCoA1", "PCoA2")) |>
  bind_cols(env) |>
  group_by(Caging, Timepoint) |>
  summarise(PCoA1=mean(PCoA1), PCoA2=mean(PCoA2))

coral_scores <- wascores(coral_pcoa$points, coral) |>
  as.data.frame() |>
  rownames_to_column() |>
  setNames(c("Taxa", "PCoA1", "PCoA2"))


a <- ggplot() + theme_bw() +
  geom_point(data = coral_pcoa_scores, aes(PCoA1, PCoA2, fill=Caging), shape=21) +
  geom_text(data = coral_pcoa_scores, aes(PCoA1-0.01, PCoA2+0.01, label=Timepoint/7), size=2.5) +
  geom_path(data = coral_pcoa_scores, aes(PCoA1, PCoA2, group=Caging),
            colour = "grey70", linewidth = 0.4, alpha = 0.85) +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.key = element_blank()) +
  guides(fill = guide_legend(title = NULL))

b <- ggplot() + theme_bw() +
  geom_segment(data = coral_scores, aes(x = 0, y = 0, xend = PCoA1, yend = PCoA2),
      arrow = arrow( length = unit(0.18, "cm"), type = "closed"),
      colour = "grey75", linewidth = 0.5, inherit.aes = FALSE) +
  geom_text_repel(data = coral_scores, aes(x = PCoA1, y = PCoA2, label=Taxa), size=2)

a|b
```

![](3-random-effects_files/figure-html/data-1.png)

There are 5 Timepoints, and 18 tiles in total (n=6 per treatment).

The correct structure to explore the effect of caging over time is
`Caging * Timepoint + (1 | Tile)`, where Timepoint as a fixed
within-tile factor whose effect (and its interaction with Caging) we
want to test, while Tile is a random intercept absorbing tile-to-tile
variation.

Using `distmer` with tile as a random effect, both `Timepoint` and
`Caging` treatment have a significant effect on benthic community
structure:

``` r
fm1 <- permanova(coral_bray ~ Caging + Timepoint + (1 | Tile), env, type="III", permutations=9999, seed=101)
fm1$table |> tab()
```

| Source | df | SS | MS | Pseudo-F | P(perm) | P(MC) | Unique perms |
|:---|---:|---:|---:|---:|---:|---:|---:|
| Caging | 2 | 1.871 | 0.935 | 4.894 | 0.005 | 0.006 | 9999 |
| Timepoint | 1 | 0.450 | 0.450 | 5.907 | 0.004 | 0.005 | 9999 |
| Tile (between-group stratum) | 15 | 2.867 | 0.191 | NA | NA | NA | NA |
| Residual | 71 | 5.412 | 0.076 | NA | NA | NA | NA |
| Total | 89 | 10.601 | NA | NA | NA | NA | NA |

With (1\|Tile) as a random effect, pseudo-F (4.89) is (correctly) tested
against MS 0.191 (tile-level), not against MS 0.076 (residual) as under
[`adonis2()`](https://vegandevs.github.io/vegan/reference/adonis.html).
Testing it against the residual would have inflated the apparent
significance, so the correct stratification (1\|tile) is more
conservative.

Exploring the pairwise contrasts, the `Caging` treatment operates
exactly as expected: Open and Partial treatments show no difference,
while Caged treatments differ from both Open and Partial:

``` r
permanova_contrasts(fm1, type = "pairwise", factor="Caging") |> tab()
```

| contrast | pseudo_F | df_num | df_den | denominator | p_perm | unique_perms | estimable | p_adj |
|:---|---:|---:|---:|:---|---:|---:|:---|---:|
| Caged - Open | 6.388 | 1 | 10 | MS\[between-group stratum\] | 0.008 | 9999 | TRUE | 0.024 |
| Caged - Partial | 5.045 | 1 | 10 | MS\[between-group stratum\] | 0.018 | 9999 | TRUE | 0.035 |
| Open - Partial | 0.510 | 1 | 10 | MS\[between-group stratum\] | 0.629 | 9999 | TRUE | 0.629 |

### Nested hierarchies with unequal effort

Ecological surveys are rarely balanced. `(1 | region/site/quadrat)`
states the containment; `distmer` then tests each fixed term at the
level where it actually varies, and does not require equal replication
at any level.

The survey below has eight regions, two sites in each, and either three
or two quadrats per site depending on which survey round the region fell
in. Reef status is a property of the region; wave exposure varies
between replicates inside a quadrat.

``` r
set.seed(2024)
quadrats_per_site <- c(3, 3, 2, 2, 3, 3, 2, 2)

survey <- do.call(rbind, lapply(1:8, function(region) {
  do.call(rbind, lapply(1:2, function(site) {
    do.call(rbind, lapply(seq_len(quadrats_per_site[region]), function(quadrat) {
      data.frame(region = paste0("R", region), site = paste0("S", site),
                 quadrat = paste0("Q", quadrat), replicate = 1:3)
    }))
  }))
}))

for (v in c("region", "site", "quadrat")) survey[[v]] <- factor(survey[[v]])
survey$status <- factor(
  ifelse(as.integer(survey$region) <= 4, "reference", "impacted"),
  levels = c("reference", "impacted")
)
survey$exposure <- factor(rep_len(c("sheltered", "exposed"), nrow(survey)))

region_effect <- rnorm(8, 0, 0.6)
latent <- 1.1 * (survey$status == "impacted") +
  0.7 * (survey$exposure == "exposed") +
  region_effect[as.integer(survey$region)] + rnorm(nrow(survey), 0, 0.5)
species <- sapply(1:12, function(k) {
  rpois(nrow(survey), lambda = exp(1.2 + 0.45 * k / 12 * latent))
})
colnames(species) <- paste0("sp", 1:12)
survey_bray <- vegdist(species, method = "bray")

fm_nested <- permanova(
  survey_bray ~ status + exposure + (1 | region / site / quadrat),
  survey, permutations = 999, seed = 11
)
fm_nested$table |> tab()
```

| Source         |  df |    SS |    MS | Pseudo-F | P(perm) | P(MC) | Unique perms |
|:---------------|----:|------:|------:|---------:|--------:|------:|-------------:|
| status         |   1 | 0.053 | 0.053 |   14.419 |   0.056 |    NA |          999 |
| exposure       |   1 | 0.291 | 0.291 |    7.214 |   0.001 |    NA |          999 |
| region stratum |   6 | 0.022 | 0.004 |       NA |      NA |    NA |           NA |
| Residual       |  79 | 3.186 | 0.040 |       NA |      NA |    NA |           NA |

Each term found its own stratum. `status` is constant within a region,
so it is tested against variation among regions (6 df) with whole
regions permuted. `exposure` varies inside a quadrat, so it is tested
against the residual (79 df) with observations permuted within quadrats.
Testing `status` against the residual would treat 120 replicates as 120
independent units and badly overstate the evidence.

Because effort is unequal, the question of *what* a region centroid
means matters. `distmer` averages upwards in stages: replicates onto
quadrat centroids, quadrat centroids onto their site with equal weight
per quadrat, and site centroids onto their region with equal weight per
site. A region surveyed with three quadrats per site therefore does not
count for more than one surveyed with two. Permutation respects the same
structure: a region moves only among regions whose whole sub-tree has
the same shape, so the two survey rounds are never exchanged with one
another.

Any number of fixed terms is allowed, including interactions and
continuous covariates, and `type` selects the sums-of-squares partition
among the terms that share a level.

`distinct_statistics` and `minimum_attainable_p` are recorded on every
test. The level the warning is written against is the
`distmer.resolution_alpha` option, `0.05` by default; see
[`` ?`distmer-options` ``](https://marine-ecologist.github.io/distmer/reference/distmer-options.md)
for that and the other reporting options, none of which change an
inferential result.

### Random slopes:

For a more complex example:

The
[`vegan::pyrifos`](https://vegandevs.github.io/vegan/reference/pyrifos.html)
data contain aquatic-invertebrate communities from 12 ditches sampled
repeatedly before and after insecticide application. Dose was assigned
to ditches, so ditch is the experimental unit. Here four balanced
post-treatment dates are used per ditch, where each trajectory is an
individual ditch according to dose:

``` r

data(pyrifos, package = "vegan")
ditch <- gl(12, 1, length = nrow(pyrifos))
week <- rep(c(-4, -1, 0.1, 1, 2, 4, 8, 12, 15, 19, 24), each = 12)
dose_by_ditch <- c(0.1, 0, 0, 0.9, 0, 44, 6, 0.1, 44, 0.9, 0, 6)
dose <- factor(rep(dose_by_ditch, 11))
keep <- week %in% c(1, 2, 4, 8)

community <- pyrifos[keep, ]
pyrifos_env <- data.frame(
  ditch = droplevels(ditch[keep]),
  dose = droplevels(dose[keep]),
  week = week[keep]
)
pyrifos_bray <- vegdist(community, method = "bray")

ditch_ord <- cmdscale(pyrifos_bray, k = 2)
ditch_ord <- cbind(
  pyrifos_env,
  data.frame(Axis1 = ditch_ord[, 1], Axis2 = ditch_ord[, 2])
)

ggplot() +
  geom_path(data = ditch_ord, aes(Axis1, Axis2, group = ditch, colour = dose, order = week), alpha = 0.65, linewidth = 0.6) +
  geom_point(data = ditch_ord, aes(Axis1, Axis2, group = ditch, fill = dose, order = week), shape=21, alpha = 0.85) +
  coord_fixed() +
  theme_bw() +
  scale_fill_brewer(palette="Set1") +
  scale_color_brewer(palette="Set1") +
  scale_size_continuous(range = c(1.8, 4), breaks = c(1, 4, 8)) +
  labs(x = "PCoA axis 1", y = "PCoA axis 2", colour = "Dose",
       size = "Week")
#> Warning in geom_path(data = ditch_ord, aes(Axis1, Axis2, group = ditch, :
#> Ignoring unknown aesthetics: order
#> Warning in geom_point(data = ditch_ord, aes(Axis1, Axis2, group = ditch, :
#> Ignoring unknown aesthetics: order
#> Ignoring unknown labels:
#> • size : "Week"
```

![Post-treatment Pyrifos samples connected through time within each
ditch](3-random-effects_files/figure-html/repeated-ordination-1.png)

For a balanced two-condition comparison, an early (1-2 weeks) and late
(\>3 week) period is defined for each ditch (balanced with two
observations in each ditch-period cell).

Visualising the mean PCoA axis 1 within ditch and period indicates
ditches need not share the same early-to-late change:

``` r
period_env <- transform(
  pyrifos_env,
  period = factor(
    ifelse(week %in% c(1, 2), "early", "late"),
    levels = c("early", "late")
  )
)

period_ord <- transform(ditch_ord, period = period_env$period)
period_centres <- aggregate(
  cbind(Axis1, Axis2) ~ ditch + period,
  data = period_ord, FUN = mean
)

ggplot() +
  geom_line(data = period_centres, aes(period, Axis1, group = ditch), colour = "grey65", linewidth = 0.7) +
  geom_point(data = period_centres, aes(period, Axis1, group = ditch, fill = ditch), shape=21, size = 2.4, show.legend = FALSE) +
  theme_bw() +
  scale_fill_brewer(palette="Set1") +
  scale_color_brewer(palette="Set1") +
  labs(x = "Period", y = "Mean PCoA axis 1 within ditch and period")
#> Warning in RColorBrewer::brewer.pal(n, pal): n too large, allowed maximum for palette Set1 is 9
#> Returning the palette you asked for with that many colors
```

![Early to late changes on PCoA axis 1 for each Pyrifos
ditch](3-random-effects_files/figure-html/period-plot-1.png)

Using a random intercept allows ditches to have different average
communities, while a random slope allows their within-ditch changes to
differ.

``` r
period_fit <- permanova(
  pyrifos_bray ~ period + (1 + period | ditch),
  period_env, permutations = "all", seed = 26
)

summary(period_fit) |> tab()
#> Distance-based PERMANOVA
#> Engine: random slope
#> 
#> Factors
#>  Name           Type          Levels
#>  period         Fixed          2    
#>  1+period|ditch Random, slope 12    
#> 
#> PERMANOVA table of results
#>  Source                       df SS    MS    Pseudo-F P(perm) Unique perms
#>  period                       1  0.424 0.424 6.7905   0.0005  4096        
#>  Between-subject random slope 11 0.687 0.062                              
#> 
#> Denominators
#>  Source Denominator                                 
#>  period MS[between-subject random-slope differences]
#> 
#> Estimates of components of variation
#>                  random_intercept random_slope
#> random_intercept          0.05305      0.00388
#> random_slope              0.00388      0.01552
```

| Source | df | SS | MS | Pseudo-F | P(perm) | P(MC) | Unique perms |
|:---|---:|---:|---:|---:|---:|---:|---:|
| period | 1 | 0.424 | 0.424 | 6.791 | 0 | NA | 4096 |
| Between-subject random slope | 11 | 0.687 | 0.062 | NA | NA | NA | NA |

A numeric gradient can also vary by unit. With four post-treatment
weeks, `(1 + week | ditch)` estimates one within-ditch community slope
and tests the equally weighted average slope against the variation among
ditch slopes.

``` r
continuous_fit <- permanova(
  pyrifos_bray ~ week + (1 + week | ditch),
  pyrifos_env, permutations = 999, seed = 26
)

summary(continuous_fit) |> tab()
#> Distance-based PERMANOVA
#> Engine: continuous random slope
#> 
#> Factors
#>  Name         Type          Levels
#>  week         Fixed         NA    
#>  1+week|ditch Random, slope 12    
#> 
#> PERMANOVA table of results
#>  Source                     df SS    MS    Pseudo-F P(perm) Unique perms
#>  week                       1  0.016 0.016 6.8424   0.0010  999         
#>  Between-group random slope 11 0.025 0.002                              
#>  Residual                   24 1.067 0.044                              
#> 
#> Denominators
#>  Source Denominator                    
#>  week   MS[between-group random slopes]
#> 
#> Estimates of components of variation
#>  Source   Estimate     Sq.root   
#>  slope    0.0007352657 0.02711578
#>  residual 0.0444687561 0.21087616
```

| Source                     |  df |    SS |    MS | Pseudo-F | P(perm) | P(MC) | Unique perms |
|:---------------------------|----:|------:|------:|---------:|--------:|------:|-------------:|
| week                       |   1 | 0.016 | 0.016 |    6.842 |   0.001 |    NA |          999 |
| Between-group random slope |  11 | 0.025 | 0.002 |       NA |      NA |    NA |           NA |
| Residual                   |  24 | 1.067 | 0.044 |       NA |      NA |    NA |           NA |

### Only random effects

`distmer` tests fixed effects against random strata, so a model with no
fixed term has nothing to test.
[`permanova()`](https://marine-ecologist.github.io/distmer/reference/permanova.md)
refuses it rather than returning a sums of squares table with no
statistic, which would read like a result:

``` r
permanova(pyrifos_bray ~ (1 | ditch), pyrifos_env)
#> Error:
#> ! No fixed term in the model. distmer tests fixed effects against random strata; it does not report a pure random-effects decomposition through `permanova()`. Either add the fixed factor of interest, for example `d ~ treatment + (1 | ditch)`, or call `estimate_variance_components()` on the compiled design for the variance components themselves.
```

The variance components themselves *are* available — they simply are not
a test. Call the estimator directly on the compiled design:

``` r
random_only <- compile_permanova_design(
  ~ (1 | ditch), pyrifos_env, pyrifos_bray
)
estimate_variance_components(random_only)$unrestricted
#>   residual      ditch 
#> 0.06216244 0.04924920
```

For a `/`-generated hierarchy use
[`estimate_nested_variance_components()`](https://marine-ecologist.github.io/distmer/reference/estimate_nested_variance_components.md),
which returns one component per level plus the residual and reproduces
the classical nested analysis-of-variance decomposition exactly.

### What distmer covers

`distmer` covers designs built on a single random structure — one
grouping factor, one containment chain of any depth, or one declared
treatment-by-block interaction — and, through
[`permanova_crossed_incomplete()`](https://marine-ecologist.github.io/distmer/reference/permanova_crossed_incomplete.md),
designs whose random terms are crossed rather than chained. Within
either scope, imbalance is allowed throughout, and cells may be missing
entirely.

| Formula | Design | Tested against |
|:---|:---|:---|
| (1 \| site) | repeated observations within sites | among-site variation |
| (1 \| site) with several fixed factors | split-plot / randomized block, with or without covariates | unit stratum or within-unit residual, per term |
| (1 \| region/site), (1 \| region/site/quadrat), deeper | nested hierarchy of any depth | the stratum where each term varies |
| (1 \| block) + (1 \| treatment:block) | treatment randomized within blocks | declared treatment-by-block variation |
| (1 \| unit) + (1 \| unit:subplot) | asymmetric nested design, three error strata | unit, unit-by-subplot, and residual strata |
| (1 \| region/site) + (1 \| treatment:region) + (1 \| treatment:region:site) | treatment crossed with a nested hierarchy | treatment-by-region interaction |
| (1 \| a) + (1 \| b), with or without (1 \| a:b) | two crossed random factors | the stratum carrying the matching components |
| (1 \| t:unit) + (1 \| t:unit:sub) + (1 \| t:unit:time) | nested hierarchy crossed with a repeated measure | the stratum each term varies over, by expected mean squares |
| numeric covariates alongside any of the above | covariate-adjusted designs | the residual, at one degree of freedom per slope |
| no random term at all, with empty cells | incomplete fixed factorial | the residual |
| (1 + condition \| subject) | complete two-condition repeated measures | among-subject variation in the condition difference |
| (1 + x \| group) | random slopes on a numeric covariate | among-group variation in slopes |

Unequal replication is supported at every level of a nesting chain, in
single-intercept designs, and in incomplete blocks. Empty cells reduce a
term to the dimension it can actually estimate rather than being
refused, and `nominal_df` records what a complete design would have
given. The two random-slope forms are the exception: both require
complete, equally replicated units.

Crossed random structures are handled by
[`permanova_crossed_incomplete()`](https://marine-ecologist.github.io/distmer/reference/permanova_crossed_incomplete.md),
which derives each denominator from the design’s expected mean squares
under the restricted parameterisation, and reports a Freedman-Lane
permutation p-value. `hypothesis` selects Type I, II or III sums of
squares; Type II is the default because it stays well posed when cells
are missing, where Type III can leave a main effect spanned by a term
containing it and report no estimable degrees of freedom.

### What distmer doesn’t cover (yet)

The remaining boundary is a random slope alongside anything else, and
hierarchies that branch.

| Formula | Status | Instead |
|:---|:---|:---|
| (1 + x \| g) with any other random term | refused | drop to a single random structure |
| (1 \| a/b) + (1 \| c) | refused: a chain crossed with an unrelated factor | declare the validated hybrid form, or drop the crossed term |
| (1 \| a/b) + (1 \| a/c) | refused: branching hierarchies | one containment chain only |
| ~ (1 \| g) with no fixed term | refused: there is nothing to test | estimate_variance_components(), or estimate_nested_variance_components() |
| a factor that is both fixed and random | refused | usually a misdeclared ‘/’; see the error text |
| (0 + x \| g), (x \| g), (1 + x \|\| g) | not accepted | write (1 + x \| g); uncorrelated slopes are unsupported |

Correlation structures such as AR(1) or spatial covariance, and
group-varying residual variance, are outside the formula grammar
entirely.

Refusals are deliberate. A formula that parses is not automatically a
defensible permutation design, and `distmer` declines a design rather
than assign it a denominator it cannot justify.
[`permanova_capabilities()`](https://marine-ecologist.github.io/distmer/reference/permanova_capabilities.md)
returns the same contract in machine-readable form.

## Implementation

For exact detail, permanova() performs the following sequence:

1.  It validates the dissimilarity and computes `G = -0.5 J D^2 J`. The
    same centring applies to every distance matrix; Euclidean distances
    give a positive semi-definite `G`, while signed axes from
    non-Euclidean distances remain in the analysis.
2.  It converts fixed terms into estimable model subspaces and records
    whether each random term is a main grouping factor, a nesting
    declaration, an explicit interaction, or a random slope.
3.  It checks balance, containment, crossing, cell completeness,
    confounding, and the level at which each predictor varies.
4.  For supported balanced designs it derives the exact
    expected-mean-square (EMS) denominator. For supported unbalanced
    designs it works with equally weighted unit centroids or profiles
    and never invents an unbalanced EMS.
5.  It constructs the permitted residual rearrangements, intact-unit
    permutations, or whole-unit sign flips, then reports the estimand,
    denominator, exchangeable unit, attainable permutations, and
    assumptions.

If no validated method matches the declared design, the model stops. A
formula that parses is not automatically a defensible permutation
design.

Most users only need
[`permanova()`](https://marine-ecologist.github.io/distmer/reference/permanova.md).
The functions below are useful for inspecting, validating, extending a
particular design, or understanding how `distmer` works:

| Formula pattern | Specialised implementation | Denominator or permutation |
|:---|:---|:---|
| d ~ a | permanova_unbalanced(); permanova_oneway() is the direct one-factor reference | Within-group residual; unrestricted observations |
| d ~ a + b (unbalanced, type = ) | permanova_unbalanced() | Full-model residual; term-specific Freedman-Lane |
| d ~ a + (1 \| site) | permanova_mixed() | Exact whole-plot EMS |
| d ~ a + x + (1 \| site) | permanova_mixed_covariate() | Split-plot or randomized-block covariate strata |
| d ~ a + (1 \| region/site) | permanova_nested() | Level-specific EMS; descendants move intact |
| d ~ a + x + (1 \| region/site) | permanova_nested_covariate() | Covariate-adjusted nested strata |
| d ~ a + (1 \| block) + (1 \| a:block) | permanova_crossed() | Exact a:block EMS; intact cells within block |
| d ~ a + (1 \| region/site) + (1 \| a:region) + … | permanova_hybrid() | Exact a:region EMS; intact treatment profiles |
| d ~ x + (1 + x \| group) | permanova_random_slope_continuous() | Between-group slope; whole-group sign flips |
| unbalanced random / nested / block | permanova_unequal_random_intercept() / permanova_unbalanced\_\*() | Unit-centroid inference; size- or topology-restricted permutations |

[`permanova_capabilities()`](https://marine-ecologist.github.io/distmer/reference/permanova_capabilities.md)
returns the complete machine-readable support table, including
estimands, error-control claims, and explicit exclusions.
