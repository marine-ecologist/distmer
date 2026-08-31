# distmer

### distmer: distance-based mixed-effects models with permutational inference

`distmer` extends the PERMANOVA (**Perm**utational **M**ultivariate
**An**alysis **O**f **Va**riance, [Anderson
2001](https://ecoevol.ufg.br/adrimelo/div/Anderson-2001-AustEcol_non-parametric_manova.pdf))
approach of
[`vegan::adonis2()`](https://vegandevs.github.io/vegan/reference/adonis.html)
to include hierarchical mixed-effects structures common in ecology.
`distmer` uses familiar `lmer`-style formulas for fixed effects,
experimental units, nesting, blocks, repeated observations, and random
slopes:

``` r

distance ~ treatment + (1 | site)
```

where the response (`distance`) is a dissimilarity matrix from
[`vegdist()`](https://vegandevs.github.io/vegan/reference/vegdist.html)
rather than a single variable. `distmer` applies a common projection and
mixed-design framework to the Gower-centred matrix induced by the user’s
chosen dissimilarity (see [Anderson
2001](https://ecoevol.ufg.br/adrimelo/div/Anderson-2001-AustEcol_non-parametric_manova.pdf)).
With no random effects, `distmer` reduces to ordinary PERMANOVA
identical to
[`vegan::adonis2()`](https://vegandevs.github.io/vegan/reference/adonis.html).
With a univariate Euclidean response, supported models reduce to their
classical [`stats::anova()`](https://rdrr.io/r/stats/anova.html) or
mixed-model counterparts
(e.g. [`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html))

### Installation

Install the development version from GitHub:

``` r

remotes::install_github("marine-ecologist/distmer")
```

`distmer` is currently on the `0.9.4` pre-release line leading to a
stable `1.0.0` release.

### Usage:

#### i) fixed effects

The `distmer` package works directly with distances produced by
[`vegan`](https://vegandevs.github.io/vegan/). The below example uses
the `pyrifos` data ([Van den Brink and Ter Braak
1999](https://onlinelibrary.wiley.com/doi/full/10.1002/etc.5620180207))
to test whether an insecticide (chlorpyrifos) treatment altered aquatic
invertebrate community composition:

``` r

library(distmer)
library(vegan)

data(pyrifos)

ditch <- gl(12, 1, length = 132)
week  <- gl(11, 12, labels = c(-4, -1, 0.1, 1, 2, 4, 8, 12, 15, 19, 24))
dose  <- factor(rep(c(0.1, 0, 0, 0.9, 0, 44, 6, 0.1, 44, 0.9, 0, 6), 11))

pyr.env <- data.frame(ditch, week, dose)
pyr.bray <- vegdist(pyrifos)

fit <- permanova(
  pyr.bray ~ dose,
  data = pyr.env,
  permutations = 999,
  seed = 101
)

summary(fit)
```

The result reports the ecological effect being tested, its pseudo-(F),
permutation p-value, attainable permutation space, and the
exchangeability assumptions used for inference.

#### ii) random effects

For a repeated or grouped field design, add the independent sampling
unit in the same way as a mixed-effects model:

``` r

fit2 <- permanova(
  pyr.bray ~ dose + (1 | ditch),
  data = pyr.env,
  permutations = 999,
  seed = 101
)

summary(fit2)
```

#### iii) pairwise contrasts

Pairwise contrasts (treatments and custom contrasts) can also be tested:

``` r

permanova_contrasts(fit1, type = "pairwise")
```

#### iv) dispersion

A location result should be considered alongside multivariate
dispersion. The dispersion companion uses the same formula-and-data
style:

``` r

permanova_dispersion(
   pyr.bray ~ ditch ,
  data = pyr.env,
  permutations = 9999,
  seed = 101
)
```

### Available formula:

| Ecological question | Formula |
|----|----|
| Do communities differ among groups? | `d ~ habitat` |
| Does treatment matter after an environmental gradient? | `d ~ treatment + depth` |
| What is the treatment effect with repeated samples per site? | `d ~ treatment + (1 | site)` |
| At which level of a nested survey does an effect occur? | `d ~ treatment + (1 | region/site)` |
| Was treatment randomized within blocks? | `d ~ treatment + (1 | block) + (1 | treatment:block)` |
| Do sites differ in their temporal response? | `d ~ time + (1 + time | site)` |

[`permanova()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova.md)
checks whether the declared design falls inside a validated inference
envelope. Unsupported or confounded structures stop with an explanation
instead of receiving a convenient but unjustified denominator.

The package also supports planned and pairwise contrasts, multiplicity
adjustment, Type I/II/III fixed-effect hypotheses, unbalanced sampling,
nested and crossed designs, and validated random-slope models. The
guides introduce these in ecological language before exposing
implementation details.

### vignettes

1.  [`vignette("getting-started", package = "distmer")`](https://marine-ecologist.github.io/distmer/articles/getting-started.html)
    — the distance-based mixed-model idea and a first analysis.
2.  [`vignette("random-effects", package = "distmer")`](https://marine-ecologist.github.io/distmer/articles/random-effects.html)
    — repeated sites, nested sampling, blocks, hybrids, and random
    slopes.
3.  [`vignette("dispersion", package = "distmer")`](https://marine-ecologist.github.io/distmer/articles/dispersion.html)
    — distance geometry and multivariate dispersion.
4.  [`vignette("contrasts", package = "distmer")`](https://marine-ecologist.github.io/distmer/articles/contrasts.html)
    — planned comparisons, pairwise tests, and multiplicity control.
5.  [`vignette("fixed-effects", package = "distmer")`](https://marine-ecologist.github.io/distmer/articles/fixed-effects.html)
    — Type I, II, and III fixed-effect questions.
6.  `vignette("atypical-distance-matrices", package = "distmer")` —
    Atypical distance matrices.
7.  [`vignette("validation", package = "distmer")`](https://marine-ecologist.github.io/distmer/articles/validation.html)
    — reference agreement, simulations, and reproducible reporting.

### Checks:

Supported analyses are tested against `vegan`, classical scalar
analogues, direct matrix calculations, exact permutation enumeration,
and adversarial simulation. Distance diagnostics retain and report
non-Euclidean geometry rather than silently changing it.

For the detailed evidence and current support boundary, see
[`inst/VALIDATION.md`](https://marine-ecologist.github.io/distmer/dev/inst/VALIDATION.md)
and
[`permanova_capabilities()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_capabilities.md).
