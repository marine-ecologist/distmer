# Validated support matrix

This matrix is a compact map, not a replacement for the design-specific
contracts in the vignettes. `permanova_capabilities()` exposes the same map as
a data frame. For an individual compiled design, runtime validators remain authoritative
and reject unsupported or confounded structures before inference.

| Primary implementation | Validated scale and denominator | Permutation/error contract |
|---|---|---|
| `permanova_oneway` | One-factor centroids over within-group residual | Unrestricted raw observations; single-hypothesis p-value |
| `permanova_fixed` | Sequential fixed terms over full-model residual | Term-specific residual permutations |
| `permanova_mixed` | Balanced whole-plot centroids over whole-plot residual EMS | Whole-plot residual permutations |
| `permanova_mixed_covariate` | Balanced whole-plot, subplot, interaction, and supported covariate hypotheses over their identified unit or within-unit stratum | Whole-unit or within-unit Freedman--Lane residual permutation according to the estimand |
| `permanova_unbalanced_splitplot` | Equal-weight unit-centroid whole-plot hypotheses and observation-scale within-unit hypotheses | Whole units move within exact-size strata; observations move within units; hypotheses not spanned by the restricted action reject |
| `permanova_asymmetric_mixed` | Whole-plot terms over the unit stratum and subplot terms over the unit-by-subplot stratum | Whole units move within identical subplot-count topologies; observations move within units |
| `permanova_crossed` | Declared crossed fixed terms over exact or EMS-matched strata | Canonical within-block permutation where validated; other p-values require explicit parametric opt-in |
| `permanova_crossed_incomplete` | Crossed fixed and random terms where cells may be empty or replication unequal; marginal SS over EMS-matched strata, Satterthwaite where no single stratum matches | Freedman-Lane permutation of reduced-model residuals |
| `permanova_nested` | Fixed term wholly in one balanced hierarchy stratum | Intact descendants move within the matching parent level |
| `permanova_dispersion` | Mean distance to a spatial median or centroid over a design-selected denominator | Full-model residual permutation under the location design's restrictions |
| `permanova_hybrid` | Treatment profiles over treatment-by-region interaction | Intact profiles move within regions |
| `permanova_unbalanced_hybrid` | Equal-site-weight observed treatment-by-region profiles after additive region adjustment | Intact profiles move within region-by-exact-site-topology strata; no unbalanced EMS claim |
| `permanova_unbalanced` | Explicit Type I, II, or III fixed hypothesis over full-model residual | Term-specific unrestricted Freedman--Lane residuals |
| `permanova_unequal_random_intercept` | Equally weighted cluster centroids over unit residual | Whole clusters move within exact-size strata |
| `permanova_unbalanced_nested` | Equally weighted top-unit centroids, each averaging its child centroids with equal weight, over the top-unit residual | Intact top units move within complete topology strata |
| `permanova_unbalanced_nested_multi` | Any number of fixed terms in an unbalanced `/`-generated hierarchy of any depth; each term is tested on equally weighted cascaded centroids at its own level over that level's residual, with Type I/II/III partitioning among terms sharing a level, and terms varying inside the deepest unit tested over the within-unit residual | Intact units move within parent-by-identical-subtree-topology strata, observations move within the deepest units; contrasts not spanned by the restricted action reject |
| `permanova_unbalanced_nested_lower` | Equally weighted child centroids adjusted for parent, or observations adjusted for child | Child centroids move within parent-by-size strata; observations move within child |
| `permanova_nested_covariate` | Covariate-adjusted fixed term wholly in one balanced hierarchy stratum; the covariate is charged to and tested over the observation residual | Term-specific Freedman--Lane permutation at the assigned hierarchy level with the covariate held in the reduced model |
| `permanova_unbalanced_block` | Equally weighted observed cells after block adjustment | Cell centroids move within block-by-size strata |
| `permanova_random_slope` | Population mean complete-subject difference over between-subject random-slope variation | Whole-subject sign flips under an explicit central-symmetry assumption |
| `permanova_random_slope_continuous` | Equally weighted within-group slopes of an ordinary numeric covariate over between-group random-slope variation | Whole-group sign flips of within-group slope vectors under an explicit central-symmetry assumption |
| `permanova_contrasts` | Fit-based pairwise or planned contrasts; isolated two-level contrasts retain the complete fitted design, while pooled multi-level statistics remain descriptive | Isolated tests refit the complete formula on selected levels and inherit its denominator and permutation action; pooled p-values remain `NA` |
| `permanova_marginal_contrasts` | Equal- or frequency-weighted categorical marginal contrast | Contrast-specific residual permutations and explicit multiplicity adjustment |
| `permanova_adjusted_contrasts` | Fixed contrast at explicitly supplied numeric covariate values in a numerically stable full-rank model | Contrast-specific residual permutations and explicit multiplicity adjustment; aliased or ill-conditioned fixed models reject |
| `permanova_mixed_contrasts` | Planned or pairwise experimental-unit contrast over its design-selected denominator | Contrast-specific unit residual permutations within declared strata |
| `permanova_maxT` | Fixed one-factor contrast family over within-group residual | Shared raw observations; weak FWER under the complete equal-centroid null |
| `permanova_pairwise_strong` | All pairwise complete-distribution equalities for one fixed factor over within-group residual | Intersection-specific label reassignment within connected equality components; closed-test strong FWER |
| `permanova_mixed_maxT` | Mixed experimental-unit contrast family over its design-selected denominator | Shared raw units within declared strata; weak FWER under the complete treatment null |

Pairwise companion functions construct all level pairs for the corresponding
contrast implementation. Strong subset-null control is not claimed by either
shared max-F implementation. Strong control is confined to
`permanova_pairwise_strong()`: independent fixed one-factor observations and
pairwise equality of complete distributions, with a separately valid
permutation action for every intersection. It is not a centroid-only claim
under heterogeneous distributions.

Across the pre-1.0 implementation, random slopes, partial unbalanced hybrid
site profiles, unsupported nested estimands, aliased hypotheses, disconnected
or confounded permutation actions, implicit covariate values, and undeclared
covariance components remain explicit failures. Non-Euclidean dissimilarities
retain empirical inference where the implementation declares a permutation
path; parametric p-values are suppressed.
