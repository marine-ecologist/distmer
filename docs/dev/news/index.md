# Changelog

## distmer 0.9.4

## distmer 0.9.3

Version 0.9.3 extends unbalanced nested inference to hierarchies of any
depth and to any number of fixed terms, adds calibration validators for
the mixed engines introduced in 0.9.2, and corrects a series of
numerical-conditioning defects found while doing so.

**Results change for two engines.**
[`permanova_unbalanced_nested()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_unbalanced_nested.md)
now averages child centroids with equal weight, and permutation p-values
are corrected wherever an absolute comparison tolerance had made them
too small; see “Multi-level unbalanced nested designs” and “Extremeness
comparison used an absolute tolerance” below.

### Crossed mixed designs with empty cells or unequal replication

- New engine
  [`permanova_crossed_incomplete()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_crossed_incomplete.md),
  reached automatically from
  [`permanova()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova.md)
  when a crossed mixed design is not a complete, equally replicated
  factorial. Previously such a design was refused outright:
  [`permanova_crossed()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_crossed.md)
  requires every cell present and every cell the same size, and the
  unbalanced fixed engine separately refuses rank deficiency, so a field
  experiment revisited on different schedules at different sites had no
  route through the package at all.
- **Degrees of freedom are estimable dimensions, not level-count
  products.** With `s` sites, `t` times and only `c` of the `s*t` cells
  occupied, the site-by-time interaction spans `c - s - t + 1`
  dimensions rather than `(s-1)(t-1)`. The table reports what the term
  can actually estimate and carries `nominal_df` alongside, so the loss
  is visible rather than silent; `$rank_deficient` names the affected
  terms.
- Sums of squares are marginal, under sum-to-zero coding. Each term is
  compared against the model holding every term that does **not**
  contain it. Dropping a term’s columns from the full model – the
  obvious construction – is the wrong comparison here: a term containing
  this one spans it, so with the site-by-time interaction still fitted,
  removing the site columns removes nothing and site reports zero
  degrees of freedom. In a complete design the higher-order contrasts
  are orthogonal to the margin and the question never arises, which is
  why the balanced engine can drop columns safely.
- Denominators are matched from the design’s expected mean squares under
  the restricted parameterisation, so a fixed term crossed with a random
  factor is tested over the stratum carrying the same random components
  and falls back to a Satterthwaite combination when imbalance leaves no
  single stratum matching. The restriction applies on both sides: a
  random stratum carries no component from a random-by-fixed interaction
  inside it, and a fixed term’s own expected mean square carries only
  the components whose fixed factors are exactly its own – `E(MS_bird)`
  holds the bird-by-time component but not site-by-bird-by-time. Without
  the second rule a random term would be matched against itself and
  return a pseudo-F of exactly one.
- Continuous covariates are carried as slopes rather than factors, so a
  covariate costs one degree of freedom instead of one per distinct
  value.
- Nesting is told apart from crossing when deciding which random
  components belong in a term’s expected mean square. A random term
  written `treatment:chamber` names the treatment only because chamber
  is nested inside it, whereas `treatment:chamber:time` is genuinely
  crossed with time. Reading the label literally put `time` over the
  residual instead of over `chamber(treatment) x time`; the test is now
  whether the fixed factor is determined by the term’s other variables.
- Interaction terms are matched by their variable set rather than their
  written order, since [`terms()`](https://rdrr.io/r/stats/terms.html)
  can return `chamber:time` as `time:chamber`.
- Random terms are optional. With none the engine covers the
  **incomplete fixed** design, which
  [`permanova()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova.md)’s
  Type III fixed path refuses outright (“aliased columns from incomplete
  cells”): every stratum collapses to the residual and each term is
  tested over it, while the rank-based construction reports what the
  empty cells leave estimable.
- `hypothesis` selects the sums of squares: `"II"` (default) adjusts
  each term for every term that does not contain it, `"I"` is
  sequential, `"III"` adjusts for every other term in the model. On a
  complete equally replicated design all three coincide. On an
  incomplete one they do not, and Type III is the one to be careful
  with: a term containing this one can span it, so a main effect can
  report zero estimable degrees of freedom. Type II is the default
  because it stays well posed however incomplete the design is.
- **The correctness anchor is agreement.** On a complete, equally
  replicated design both this engine and the validated
  [`permanova_crossed()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_crossed.md)
  apply, and they return identical degrees of freedom, sums of squares
  and pseudo-F to machine precision. Everything else is that same
  construction carried into designs the balanced engine refuses.

### Multi-level unbalanced nested designs

- `permanova(d ~ terms... + (1 | region/site/quadrat))` now works when
  replication is unequal at more than one level, via a new engine
  [`permanova_unbalanced_nested_multi()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_unbalanced_nested_multi.md).
  Previously only a single parent/child pair was supported, and deeper
  unbalanced chains were refused.
- The engine takes **any number of fixed terms**, closing an asymmetry
  in the package: multi-factor support already existed for an unbalanced
  `(1 | unit)` design through
  [`permanova_unbalanced_splitplot()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_unbalanced_splitplot.md),
  but every unbalanced nested engine accepted exactly one categorical
  factor, so `~ A + B + (1 | reg/site)` or
  `~ A * B + (1 | reg/site/quad)` was refused. Each term is now assigned
  to the shallowest hierarchy level within which it is constant and
  tested there, with Type I/II/III partitioning among the terms sharing
  a level. A term that varies inside the deepest unit – typically an
  observation-level covariate – is tested over the within-unit residual
  with residuals permuted within those units. Verified against
  independent `lm` references at every level: region, site, and
  observation, including interactions, and against
  [`car::Anova`](https://rdrr.io/pkg/car/man/Anova.html) for the Type
  III partition.
- Permutation moves intact units within their own parent and only across
  units whose complete subtree shape matches, generalising the
  parent/child engine’s sorted-child-count rule. Treatment contrasts not
  spanned by that restricted action reject rather than being tested
  invalidly.
- **Behaviour change.**
  [`permanova_unbalanced_nested()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_unbalanced_nested.md)
  now uses the same cascaded aggregation, so a parent centroid averages
  its child centroids with equal weight instead of averaging the
  parent’s observations directly. The two agreed already when children
  were equally replicated; under unequal child replication the pseudo-F,
  and any p-value derived from it, will differ from 0.9.1. The
  observation-weighted parent centroid let a child carrying more
  replicates pull its parent, which is the imbalance the equal-weight
  estimand exists to remove, and it disagreed with the deeper chains.
  Top-level unbalanced-nested contrasts inherit the same change. Engines
  with no level between the aggregation unit and the observations
  ([`permanova_unequal_random_intercept()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_unequal_random_intercept.md),
  [`permanova_unbalanced_block()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_unbalanced_block.md),
  and the child/observation-level nested engine) are unaffected, since
  for them the two conventions coincide.
- Because that change moves published numbers and a NEWS entry is easy
  to miss, an affected fit now says so in its own warnings. The note is
  raised only when children really are unequally replicated within a
  parent; when they are not, the old and new conventions coincide and
  there is nothing to report.

### Denominators derived from the design’s expected mean squares

- [`permanova_asymmetric_mixed()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_asymmetric_mixed.md)
  applied the balanced-design rule – whole-plot terms over the unit
  stratum, sub-plot terms over the unit-by-sub-plot stratum – to designs
  that are not balanced. The package elsewhere promises never to invent
  an unbalanced expected mean square; this engine did so implicitly, and
  on a reef survey with island sizes from 2 to 240 sites it put six of
  seven fixed terms out by 25 to 95 per cent.
- A new `denominator` argument (default `"ems"`) derives the expected
  mean squares of the numerator and of every candidate stratum from the
  actual design by the trace method, then constructs the ratio from
  them. Where the match would need a negative coefficient on the
  residual, the surplus is moved to the **numerator** rather than
  truncated, keeping every coefficient non-negative, and Satterthwaite
  degrees of freedom apply to both sides – so a numerator df need not be
  a whole number. `denominator = "stratum"` restores the previous rule.
- Balanced designs give an exact single-stratum match and are unchanged.
- The interaction covariance component is orthogonalised against the
  unit main effect. Without it no random component is attributable to a
  single stratum, the construction has no solution, and the engine
  silently falls back.
- Validated against a published PRIMER workspace for an unbalanced
  three-stratum design: the construction reproduces PRIMER’s pseudo-F,
  numerator df and denominator df for every fixed term to all printed
  digits, and the trace-derived EMS coefficients match PRIMER’s to five
  significant figures. Within the engine the mean absolute log ratio
  against PRIMER falls from 0.335 to 0.037.

### Contrasts within the cells of other factors

- [`permanova_contrasts_within()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_contrasts_within.md)
  contrasts the levels of one fixed factor inside every combination of
  others, as PRIMER’s “within level” pairwise tables do.
- The denominator cannot be inherited from the parent fit. Inside one
  cell the lower strata often collapse – an ecoregion sampled on two
  islands leaves the unit-by-sub-plot stratum with no degrees of freedom
  – so each pair rebuilds the declared random strata on its own subset,
  drops those with no degrees of freedom, and matches the contrast’s
  expected mean square against the survivors. Three pairs inside one
  cell can therefore carry three different denominators and three
  different Satterthwaite degrees of freedom, which no single fitted
  model can produce.
- Residuals are permuted inside the experimental unit, so observations
  move between the levels being contrasted.
- **Expected mean squares are restricted, as PRIMER’s are.** PRIMER
  states the convention in its workspace header – “Fixed effects sum to
  zero for mixed terms” – and under that restricted model the mean
  square of a random stratum carries no variance component from a
  random-by-fixed interaction nested inside it: `E(MS_island)` contains
  no island-by-depth component, though `E(MS_island x depth)` does.
  distmer previously used the unrestricted accounting and then tried to
  recover the zero by orthogonalising the interaction’s `Z` against the
  containing stratum, which is the wrong instrument: it changes that
  component’s covariance everywhere it appears, including in the
  contrast’s own expected mean square, and left the inner coefficient
  too small by up to 18 per cent. The zero now goes where it belongs, in
  the stratum’s row of the accounting, and a component’s covariance is
  `ZZ'` as declared.
- **Permutation happens at the level of the denominator’s finest
  stratum.** When a denominator is carried by a grouping stratum rather
  than the residual, the exchangeable entity is a whole unit of that
  stratum – an island-by-depth cell – not an observation. Permuting
  observations inside the island instead treats every sample in a cell
  as an independent replicate; that is pseudo-replication, and it made
  p-values too small in 19 of 20 affected contrasts. Units of unequal
  size cannot be permuted as rows, so the exchange is done on unit
  centroids. This costs nothing: any projector constant within units
  satisfies `tr(H L G L') = tr(H G)`, so the observed statistic is
  unchanged and only the null moves. Residuals are permuted under a
  reduced model holding the coarser strata, weighted by unit size.
- **The numerator’s surplus residual reaches the permutation.** Where
  the matched weights carry more than one unit of residual, the
  construction moves the surplus onto the numerator rather than allowing
  a negative weight, so the statistic is
  `(numerator MS + surplus * residual MS) / denominator`. The
  permutation was computing the numerator mean square alone, leaving the
  null on a different scale from the observed statistic it was compared
  against. Both permutation paths now carry the surplus, and the
  unit-level path asserts that the unpermuted configuration reproduces
  the reported statistic exactly, so the two cannot drift apart again
  unnoticed.
- Validated against the same PRIMER workspace, now on **all 27
  contrasts** in the complete within-level table rather than one cell.
  Every t statistic, every denominator coefficient on all three strata,
  and every Satterthwaite degree of freedom reproduces to PRIMER’s
  printed precision, and the calls at the five per cent level are
  identical throughout.

### Known gap: the within-cell null over a grouping stratum

Where the residual carries a within-cell denominator, the exchangeable
unit is the observation and the permutation p-values match PRIMER’s to
within 0.007.

Where a *grouping* stratum carries it, the exchangeable unit is a whole
unit of that stratum. distmer now permutes the right space – PRIMER’s
own unique-perm counts name it exactly at 720 = 6! for a cell with six
occupied island-by-depth units, and distmer finds the same 720 – but the
distribution over that space is not yet identical to PRIMER’s. Gaps
reach 0.157 on 18 such contrasts, mean 0.037, with distmer the more
conservative of the two in 15 of 18.

Fixing the surplus-residual defect above moved these p-values *away*
from PRIMER, not toward it: the closer agreement beforehand was an
artifact of the null being on the wrong scale. Five variants of the
statistic and the reduced model were tested against the nine contrasts
where both sides can be enumerated exactly, and none accounts for the
difference, so the cause is not one of those choices. The tests bound
this discrepancy rather than asserting agreement. -
[`permanova_contrasts()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_contrasts.md)
reaches this family directly, via a new `within =` argument. Naming
several factors in `factor` is shorthand for the same thing: the first
is contrasted inside combinations of the rest, and the reading taken is
reported. Several factors cannot mean contrasting their crossed cells –
two cells of a crossed combination leave the design’s lower strata with
no degrees of freedom – so the shorthand is unambiguous. `within` is
refused for the non-pairwise families, which need a contrast matrix over
all levels that a single cell need not carry, and for
`method = "pooled"`, which has no meaning when every pair derives its
own denominator. - The result prints grouped by cell, so a within-level
table can be read as published.

### Whole-plot terms were untestable on unequally replicated designs

- [`permanova_asymmetric_mixed()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_asymmetric_mixed.md)
  permutes whole units by swapping blocks of observation indices, which
  is only defined between units holding the same counts in the same
  sub-plot classes. Real surveys give every unit a distinct count
  vector: in a 33-island reef dataset all 33 signatures were unique, so
  no two islands were exchangeable, the only attainable arrangement was
  the identity, and **every whole-plot term returned `p = 1.0000`
  regardless of its pseudo-F**. A term with pseudo-F of 26 was reported
  as p = 1.
- A new `whole_plot` argument decides what happens. `"observation"`
  keeps the intact-unit action at the observation level, weighting units
  by how many observations they carry – the estimand the `IS(unit)`
  stratum defines. `"centroid"` tests whole-plot terms on equally
  weighted unit centroids and permutes units freely, so a p-value exists
  for any design. `"auto"` (the default) uses the observation-level
  action whenever any two units are exchangeable and falls back to
  centroids otherwise.
- The two agree **exactly** when units are equally replicated, so no
  existing balanced analysis moves; a regression test pins that
  bit-for-bit. Under imbalance they target different estimands, and a
  fit that fell back says so in the affected terms’ warnings. The action
  taken is recorded on the fit as `$whole_plot`.
- The unit-level denominator is the residual among unit centroids after
  every whole-plot term, so all whole-plot tests share one denominator
  exactly as they do at the observation level.
- Sub-plot terms are untouched: they keep the observation-level action
  and the unit-by-sub-plot stratum.

### Permutation resolution reaches the asymmetric engine

- [`permanova_asymmetric_mixed()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_asymmetric_mixed.md)
  computed the resolution diagnostic and then discarded it, so the
  `p = 1.0000` above arrived with no explanation. It now reports
  `distinct_statistics` and `minimum_attainable_p` on every test spec
  and raises the resolution warning, which names exactly this failure
  (“the permutation distribution holds only 1 distinct statistic”).
  Progress bars in this engine were also unlabelled, so a multi-term fit
  drew several anonymous bars; each is now labelled with its term and
  denominator stratum.

### Calibration validators for the newer mixed engines

- Five engines had unit tests but no entry in the validation harness.
  Each now has a `simulate_*_validation()` function pinned to an
  independently computed classical reference, and
  [`run_permanova_validation()`](https://marine-ecologist.github.io/distmer/dev/reference/run_permanova_validation.md)
  checks that agreement:
  [`permanova_mixed_covariate()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_mixed_covariate.md)
  against `aov(y ~ tr * sub + Error(unit))`;
  [`permanova_unbalanced_splitplot()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_unbalanced_splitplot.md)
  against an `lm` on equally weighted unit centroids;
  [`permanova_nested_covariate()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_nested_covariate.md)
  against `SS(factor | covariate)` over the factor-only nested stratum;
  [`permanova_asymmetric_mixed()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_asymmetric_mixed.md)
  against `aov(y ~ EC * DE + Error(ISLAND / DE))`; and
  [`permanova_random_slope_continuous()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_random_slope_continuous.md)
  against a one-sample test on per-group slopes. All agree to better
  than 1e-13.
- Type-I error at 300 simulations is within Monte-Carlo error of the
  nominal 0.05 for all five (0.037-0.060; every 95% interval covers
  0.05).

### Restricted permutation with a singleton stratum

- `.stratified_permutation_plan()` used `sample(position)`, which
  permutes `seq_len(x)` when `x` is a single number. A permutation
  stratum containing exactly one unit therefore aborted the fit (“number
  of items to replace is not a multiple of replacement length”) unless
  that unit happened to be the first observation. Singleton strata are
  common in unbalanced data – one unit with a unique replicate count, or
  one hierarchy unit with a unique topology – and every engine using
  restricted permutation was affected, including
  [`permanova_unequal_random_intercept()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_unequal_random_intercept.md),
  [`permanova_unbalanced_splitplot()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_unbalanced_splitplot.md),
  [`permanova_asymmetric_mixed()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_asymmetric_mixed.md),
  and the unbalanced nested engines. The stratum is now indexed
  explicitly and a singleton is correctly held fixed.

### Nesting chains crossed with unrelated random factors

- [`permanova()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova.md)
  refused designs such as `~ tr + (1 | a/b) + (1 | c)` only after
  routing them into the crossed engine, which reported an internal
  “Phase 3” message about `/`-generated nesting. The dispatcher now
  refuses the structure directly, names the nesting chain involved, and
  points to either the nested engine or the validated crossed-nested
  hybrid form.

### Degenerate results from a model with no testable fixed term

- A design whose fixed terms are all absorbed by the random structure
  produced `permutation$used = -Inf` and an internal warning (“no
  non-missing arguments to max”) instead of a clean zero. Three of the
  four sites taking the maximum over an empty set of test specifications
  lacked the guard the fourth already had; all four now agree.

### Extremeness comparison used an absolute tolerance

- Permutation p-values counted a statistic as extreme when it reached
  `observed - tol` with `tol = 1e-10` **absolute**. The identity
  permutation reproduces the observed statistic only up to
  floating-point rounding, and at a pseudo-F of around 1e4 that rounding
  is roughly 1e-9 – larger than the tolerance. The permutation that
  should always tie with the observed value then failed to count, and an
  exhaustive test could report **p = 0**, which no permutation test can
  legitimately produce. In one four-unit design the reported p-value was
  0 where the correct exact value is 1/3.
- Widening the slack alone was not enough. The observed statistic and
  its permuted counterparts are reached by different arithmetic paths,
  and under severe conditioning – a pseudo-F of 1e12, say – the tie
  between them is lost by more than any fixed relative slack. Under
  exact enumeration the identity permutation is a member of the set and
  always ties, so `p >= 1/N` holds as a theorem; that bound is now
  enforced directly, which removes the whole class regardless of how
  well conditioned the design is.
- The comparison also now scales with the magnitude of the observed
  statistic in every engine, using the same `sqrt(.Machine$double.eps)`
  floor the package already applies to numerical-equality checks.
  p-values are unchanged on ordinary scales – verified identical to the
  previous idiom across 200 randomised cases – and are corrected upward
  wherever the statistic was large enough for the old slack to fail.
  Designs with well-separated groups are the ones affected, so any
  previously reported p-value of exactly 0 should be recomputed.

### Numerically-zero guards are now scale free

- Sixteen guards across twelve engines refused a fit when a mean square
  was “numerically zero”, comparing it with an **absolute**
  `tol = 1e-10`. A mean square carries the squared units of the
  dissimilarity, so that test was not scale free: the same analysis was
  accepted on a response of ordinary magnitude and refused after
  multiplying the response by 1e-4. Data already on a small scale –
  proportions, or an already-normalised dissimilarity – could be
  rejected with “the … mean square is numerically zero” when nothing was
  wrong with it. This was reported as a known limitation in 0.9.3.
- Each guard now measures the mean square against the total inertia of
  the same Gower matrix, so the question becomes “is this stratum
  negligible relative to the variation in these data”, which is what the
  guard is for and does not depend on the units of the response.

### Consistent refusal of a degenerate denominator

- [`permanova_oneway()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_oneway.md)
  and
  [`permanova_fixed()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_fixed.md)
  divided by a numerically zero residual and returned an infinite
  pseudo-F, while every other engine refused. Both now refuse with the
  same message. Only the observed statistic is guarded, so a single
  degenerate permutation cannot abort an otherwise valid run.
- [`permanova_oneway()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_oneway.md)
  gains a `tol` argument for that check.

### Permutation resolution is now reported

- A restricted permutation scheme can attain many arrangements yet
  produce only a handful of distinct statistics: permuting labels inside
  a treatment group leaves the statistic untouched. Two groups of two
  attain 24 arrangements but only three distinct pseudo-F values, so the
  smallest p-value the design can return is a third – however large the
  effect. Reporting the count of attainable arrangements alone hid this.
- Test specifications now carry `distinct_statistics` and
  `minimum_attainable_p`, and a design that cannot reach the chosen
  level however large its effect says so in its warnings. The level is
  the new `distmer.resolution_alpha` option (default `0.05`); it is a
  reporting threshold only, and the two numbers are recorded whatever it
  is set to.
- A new
  `?\`distmer-options\``help page documents every option that changes how the package reports --`distmer.resolution_alpha`,`distmer.progress`, and`distmer.progress_min\`.
  None of them changes an inferential result. The progress options
  existed since 0.9.2 but were only described in an internal comment.
- The p-value computation itself, previously duplicated at two dozen
  sites, is now a single shared helper, which is what let both of the
  above be fixed in one place.

### Models with no fixed term

- [`permanova()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova.md)
  now refuses a formula whose only terms are random, instead of
  returning a sums-of-squares table with no statistic in it. The refusal
  names the estimator that does decompose a pure random structure –
  [`estimate_nested_variance_components()`](https://marine-ecologist.github.io/distmer/dev/reference/estimate_nested_variance_components.md)
  for a `/`-generated chain,
  [`estimate_variance_components()`](https://marine-ecologist.github.io/distmer/dev/reference/estimate_variance_components.md)
  otherwise – and suggests a concrete repair using the user’s own
  grouping factor.
- [`estimate_nested_variance_components()`](https://marine-ecologist.github.io/distmer/dev/reference/estimate_nested_variance_components.md)
  no longer requires a fixed term. The requirement came from a validator
  shared with the testing path, where it is correct; the estimator
  itself never used the fixed model. With no fixed terms each stratum is
  simply left unadjusted, which is the correct pure random-effects
  decomposition, and the result now reproduces the classical nested
  analysis-of-variance components exactly.
- The “Only random effects” section of the random-effects vignette
  previously claimed an error was thrown but illustrated it with
  `~ (1 + period)`, which parses as a *fixed* term and would not have
  errored. It now shows the real refusal and the variance-component
  route.

### Removal of abandoned internal helpers

- Six never-called internal functions (`.fix1_partial_term`,
  `.fix2_distinct_perms`, `.fix3_confounded`, `.fix3_no_test_row`,
  `.fix4_classify_grouping`, `.fix5_no_fixed_message`) have been
  deleted. Each duplicated a policy the live code already implements, in
  most cases less safely: the grouping classifier counted distinct
  levels per group where the live classifier uses projection ranks, and
  the no-test table row used column names that do not match any table
  the package builds.

### Capability registry completeness

- The registry returned by
  [`permanova_capabilities()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_capabilities.md)
  was missing five routable engines (`permanova_mixed_covariate`,
  `permanova_unbalanced_splitplot`, `permanova_nested_covariate`,
  `permanova_asymmetric_mixed`, and
  `permanova_random_slope_continuous`). All are now registered, along
  with the new multi-level nested engine, and `inst/SUPPORT-MATRIX.md`
  names each one.
- A release-audit test now reads the engine names out of
  [`permanova()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova.md)’s
  own dispatcher and fails if any routable engine is absent from the
  registry, so the contract cannot drift again. The registry-size checks
  in the validation harness and unit tests no longer pin an exact count.

### Documentation

- The random-effects vignette gains two sections. “Nested hierarchies
  with unequal effort” works a `(1 | region/site/quadrat)` survey in
  which quadrat effort differs between survey rounds, showing each fixed
  term being tested at the level where it varies and explaining what a
  centroid means when effort is unequal. “What a permutation design can
  and cannot resolve” uses the same fit to show a restricted scheme
  attaining 576 arrangements but only 18 distinct pseudo-F values, so
  that the p-value cannot fall below 0.056 whatever the effect size –
  and points at
  [`` ?`distmer-options` ``](https://marine-ecologist.github.io/distmer/dev/reference/distmer-options.md)
  for the reporting level.
- Neither the multi-level nested engine nor the resolution diagnostic
  had any user-facing documentation before this.
- The same vignette now opens with two scope tables: the random
  structures that are covered, with what each fixed term is tested
  against, and those that are not yet covered, with what to reach for
  instead. The boundary is two independent random structures at once –
  one hierarchy of any depth is well covered, but a slope alongside
  another random term, or two crossed hierarchies, are not. Crossed
  random factors are listed as a partial case: the pseudo-F, EMS and
  variance components are available, only the permutation p-value is
  not.

## distmer 0.9.2

Version 0.9.2 consolidates the complete post-0.9.1 development tree: the
low-rank Freedman–Lane refactor, progress reporting, asymmetric
three-stratum mixed models, balanced and unbalanced multi-factor
split-plots, and the fit-based contrasts API. The adversarial
corrections below are part of the same reviewed source snapshot.

### Adversarial inference hardening

- Exact asymmetric-mixed and split-plot requests now enumerate their
  complete restricted spaces instead of returning zero permutations and
  `NA` p-values.
- Incomplete asymmetric units move only within identical subplot-count
  topologies, and unequal split-plot units move only within exact-size
  strata.
- Fit-based isolated contrasts retain the complete original formula,
  random effects, nuisance terms, denominator, and permutation action.

### Low-rank Freedman–Lane permutation

- Every residual-permutation engine now shares one kernel
  (`.freedman_lane_kernel()`) instead of forming the permuted Gower
  matrix `A %*% G %*% t(A)` explicitly. Expanding that product leaves
  the reduced-model projector and the cross terms as fixed matrices, so
  only `R G R` is re-indexed by each permutation; factoring each
  hypothesis projector as `W diag(s) W'` then reduces a permutation to
  one narrow matrix product with `r` columns, where `r` is the
  hypothesis rank. Residual projectors, whose rank is close to `n`, are
  carried through their complement using the fact that permutation
  similarity leaves `tr(R G R)` unchanged. Cost per permutation falls
  from `O(n^3)` to `O(n^2 r)`: roughly 11x faster at `n = 160`, 20x at
  `n = 320`, and 18x at `n = 1000`. Results are unchanged to machine
  precision (relative difference below 1e-13) and the random-number
  stream is untouched, so seeded p-values reproduce exactly. Engines
  covered:
  [`permanova_mixed()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_mixed.md),
  [`permanova_nested()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_nested.md),
  [`permanova_unbalanced()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_unbalanced.md),
  [`permanova_mixed_covariate()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_mixed_covariate.md),
  [`permanova_nested_covariate()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_nested_covariate.md),
  [`permanova_unequal_random_intercept()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_unequal_random_intercept.md),
  [`permanova_unbalanced_nested()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_unbalanced_nested.md),
  [`permanova_unbalanced_nested_lower()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_unbalanced_nested_lower.md),
  [`permanova_unbalanced_block()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_unbalanced_block.md),
  [`permanova_unbalanced_splitplot()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_unbalanced_splitplot.md),
  [`permanova_asymmetric_mixed()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_asymmetric_mixed.md),
  and every contrast family built on them.
- [`permanova_dispersion()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_dispersion.md)
  evaluates its rank-one permuted response as the quadratic form `z'Hz`
  rather than building the `n`-by-`n` outer product.

### Permutation progress reporting

- Long permutation runs now draw a text progress bar labelled with the
  term or contrast being permuted. Reporting is off outside an
  interactive session and is suppressed below 200 permutations, so
  scripts, tests, vignettes, and `R CMD check` are unaffected. Set
  `options(distmer.progress = TRUE)` to force it on (or `FALSE` to
  silence it), and `options(distmer.progress_min = )` to change the
  length at which a run is worth reporting.

### Asymmetric nested mixed designs (three error strata)

- `permanova(distance ~ fixed... + (1 | unit) + (1 | unit:subplot))` now
  routes to a new engine,
  [`permanova_asymmetric_mixed()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_asymmetric_mixed.md),
  reproducing PRIMER’s asymmetric nested mixed design with three
  observation-level error strata: whole-plot fixed effects (not
  involving the sub-plot factor) tested over the unit stratum
  `IS(unit)`; sub-plot fixed effects tested over the unit-by-sub-plot
  interaction stratum `subplot x IS(unit)`; and both random strata
  tested over the residual. Fixed-term sums of squares are standard
  full-model Type I/II/III; the whole-plot vs sub-plot classification
  only selects the denominator. Degrees of freedom are derived from the
  actual (possibly incomplete) design. For the Caribbean benthos example
  this is
  `bray ~ DEPTH * ECOREGION * POP_STATUS + (1 | ISLAND) + (1 | ISLAND:DEPTH)`,
  giving `IS(EC x PS)` (whole-plot denominator) and `DE x IS(EC x PS)`
  (sub-plot denominator) – the ST5 structure. Validated to reproduce
  `aov(y ~ DE*EC*PS + Error(ISLAND/DE))` exactly for every fixed term F,
  the stratum degrees of freedom, and the residual. (The two-stratum
  [`permanova_unbalanced_splitplot()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_unbalanced_splitplot.md)
  remains for designs without a random unit-by-sub-plot interaction.)

### Unbalanced multi-factor split-plot designs

- `permanova(distance ~ f1 + f2 + ... + (1 | unit))` now works when the
  units are **unequally replicated**, via a new engine
  [`permanova_unbalanced_splitplot()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_unbalanced_splitplot.md).
  Whole-plot factors (constant within each unit) are tested on **equally
  weighted unit centroids** over the between-unit stratum – the
  unit-weighted estimand of
  [`permanova_unequal_random_intercept()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_unequal_random_intercept.md)
  generalised to any number of fixed factors and Type I/II/III sums of
  squares – with size-restricted whole-unit permutation. Within-unit
  (sub-plot) factors and between-by-within interactions are tested over
  the within-unit residual with within-unit residual permutation. This
  is the distance-based analogue of PRIMER’s asymmetric nested design
  (e.g. `IS(EC x PS)` as the whole-plot stratum and `DE x IS(EC x PS)`
  as the sub-plot residual) for unequally sampled units. The
  [`permanova()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova.md)
  front door routes a balanced design to
  [`permanova_mixed_covariate()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_mixed_covariate.md)
  and an unbalanced one here; a lone whole-plot factor still uses
  [`permanova_unequal_random_intercept()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_unequal_random_intercept.md).
  Validated to reproduce, exactly, an `lm` on equal-weight unit
  centroids for the whole-plot terms, a unit-blocked `lm` for the
  sub-plot terms, the balanced split-plot engine when balanced, and the
  single-factor unit-centroid engine.

### Multi-factor split-plot and randomized-block mixed designs (no covariate needed)

- The single-random-intercept engine no longer requires exactly one
  fixed factor or that the factor be constant within the unit.
  `permanova(distance ~ f1 + f2 + ... (+ interactions) + (1 | unit))`
  now supports any number of fixed factors, each classified relative to
  the random unit and tested over the matching stratum: **whole-plot**
  factors (constant within every unit) over the between-unit stratum;
  **within-unit** sub-plot / randomized-block factors (fully crossed
  within every unit) and **between-by-within** interactions over the
  within-unit residual. A **partially confounded** factor (neither
  constant nor fully crossed within the unit) is refused rather than
  given a wrong denominator. A continuous covariate remains optional and
  is handled as before, including a between-unit covariate (constant
  within units), which is treated as a whole-plot term.
- [`permanova_mixed_covariate()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_mixed_covariate.md)
  (still the engine, now covariate-optional) is documented as the
  general split-plot / randomized-block / ANCOVA engine.
  [`permanova()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova.md)
  routes the simple single-whole-plot-factor case to the exact-EMS
  [`permanova_mixed()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_mixed.md)
  and everything richer to this engine. Validated to reproduce
  `aov(y ~ ... + Error(unit))` and the blocked linear model exactly for
  whole-plot factors, sub-plot factors, and between-by-within
  interactions. This realises and validates an earlier internal
  split-plot design sketch, now removed.

### Fit-based contrasts API

- [`permanova_contrasts()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_contrasts.md)
  now operates on a fitted
  [`permanova()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova.md)
  object (rather than a compiled design), with contrast families via
  [`contrast_matrix()`](https://marine-ecologist.github.io/distmer/dev/reference/contrast_matrix.md)
  and two testing methods: `"isolated"` (a validated two-group subset
  test) and `"pooled"` (whose pseudo-F equals the classical
  constrained-scalar contrast F, verified to machine precision).
  `permanova_pairwise()` is removed –
  `permanova_contrasts(fit, type = "pairwise")` replaces it. The
  validation suite and capability registry were updated to the new API.
  Shared internal helpers (`%||%`, `.stopf`, `.warnf`) are consolidated
  in one `utils.R`.

## distmer 0.9.0.9000

### Versioning reset before the first stable release

- The active package version is now `0.9.0.9000`: a pre-release
  development line leading to the first stable `1.0.0` release.
- Earlier version-like numbers were internal implementation checkpoints,
  not published package releases. Their detailed notes remain below for
  audit traceability but are explicitly labelled as unreleased
  checkpoints.
- Four-component `.9000` versions are development snapshots. The public
  release will use `1.0.0`, with later releases following semantic
  versioning.

*The source previously carrying the internal `2.9.1.9000` label was
renumbered to `0.9.0.9000`; this is a version-policy correction, not a
rollback of code or validated capability.*

### Mixed-covariate and asymptotic-inference hardening

- Rebuilt the documentation as a seven-guide ecological workflow: an
  explicit distance-based mixed-model introduction, fixed effects,
  pairwise contrasts, distance and dispersion diagnostics, unbalanced
  designs, random effects, and reproducible validation. Worked examples
  now use
  [`vegan::dune`](https://vegandevs.github.io/vegan/reference/dune.html),
  [`vegan::mite`](https://vegandevs.github.io/vegan/reference/mite.html),
  and
  [`vegan::pyrifos`](https://vegandevs.github.io/vegan/reference/pyrifos.html),
  with `ggplot2` figures. Implementation routing is kept as a secondary
  reference in the random-effects guide.

- Mixed and nested covariate projectors now use scale-aware stratum
  projection, so structurally zero between/within components cannot
  acquire numerical rank. Single-random-intercept models assign a
  numeric term that is constant within groups to the between-group
  stratum. Nested covariate inference rejects a numeric term with no
  deepest-stratum variation because that higher-level term requires a
  denominator outside the supported envelope.

- Nested covariate Type I/II/III hypotheses now apply to covariate main
  effects and slope interactions as well as factors. Their Freedman–Lane
  reduced model retains every term used to adjust the observed
  numerator.

- Balanced random-intercept inference now honours Type I/II/III
  fixed-effect hypotheses instead of coupling the numerator partition to
  the random-effect denominator choice.

- `P(MC)` now derives its eigenvalue spectrum from each test’s
  orthogonal numerator-plus-denominator stratum. Variation outside that
  stratum can no longer change an upper-level mixed-model p-value.

- `permanova_mixed_covariate(permutations = "all")` now enumerates the
  actual whole-group or within-group action, respects `max_exact`, and
  reports its finite attainable space.

- Continuous random-slope sign-flip inference remains available with
  unequal covariate information. The parametric `F(1, ng - 1)` reference
  is now limited to rank-one Euclidean slopes with equal within-group
  `SSx`; the slope variance component subtracts
  `MS[Residual] * mean(1 / SSx)`.

- The source release audit scans all Rd files and the current
  seven-vignette series. Development handoff and readiness documents now
  follow the active version rather than the retired Phase 28 tree.

*Internal checkpoint label: `2.9.0.9000` (unreleased).*

### Type I / II / III sums of squares for the fixed part of mixed designs

- The `type` argument of
  [`permanova()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova.md)
  now applies to the covariate-adjusted mixed engines, not only pure
  fixed-effect models.
  [`permanova_mixed_covariate()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_mixed_covariate.md)
  and
  [`permanova_nested_covariate()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_nested_covariate.md)
  gained a `hypothesis` argument (`"I"` sequential, `"II"`
  marginal/containment-respecting, `"III"` partial), threaded from the
  front-door `type`. The sum-of-squares type changes only how each fixed
  term’s numerator reduced model is formed; the random-structure
  denominator is unchanged. Types coincide for an orthogonal balanced
  fixed structure and differ only when the fixed terms are
  non-orthogonal (crossed factors with unequal cell counts). Validated
  against [`car::Anova`](https://rdrr.io/pkg/car/man/Anova.html): for a
  block design with an unbalanced within-block `A*B` structure, distmer
  reproduces the Type I, Type II, and (through the front door’s
  automatic sum-to-zero contrasts) Type III F ratios exactly. The front
  door no longer warns that `type` is ignored for these engines.

*Internal checkpoint label: `2.8.0.9000` (unreleased).*

### Covariate adjustment for single random-intercept designs (split-plot / RCB ANCOVA)

- `permanova(distance ~ factor + covariate + (1 | group))` now routes to
  a new engine,
  [`permanova_mixed_covariate()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_mixed_covariate.md),
  extending covariate adjustment beyond the balanced nested engine to
  the single-random-intercept case. Each fixed term is classified by
  whether it varies within the grouping unit: between-unit (whole-plot)
  factors are tested over the between-group random-intercept stratum,
  while within-unit terms – sub-plot or randomized-block factors,
  observation-level covariates, and covariate-by-factor slope
  interactions – are tested over the within-group residual. This is the
  distance-based split-plot and randomized-block analysis of covariance.
  Unlike
  [`permanova_nested_covariate()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_nested_covariate.md),
  it removes the covariate’s between-unit component from the factor
  denominator and adjusts each whole-plot factor for the covariate’s
  *between-unit* slope (which differs from the pooled slope), so a
  covariate with between-unit structure is handled correctly without a
  manual
  [`decompose_covariate()`](https://marine-ecologist.github.io/distmer/dev/reference/decompose_covariate.md)
  step. Validated to reproduce the whole-plot factor test against
  `aov(y ~ covariate + factor + Error(group))` and the within-unit
  (sub-plot factor, covariate, slopes) tests against the corresponding
  blocked linear model, exactly, including degrees of freedom.

*Internal checkpoint label: `2.7.0.9000` (unreleased).*

### Continuous-covariate random slopes

- `permanova(distance ~ x + (1 + x | group))` with a numeric `x` now
  routes to a new engine,
  [`permanova_random_slope_continuous()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_random_slope_continuous.md),
  for the covariate-by- random model of Anderson, Gorley & Clarke (2008,
  Fig. 1.48b). The fixed slope has no exact single-source denominator in
  the full observation space (its expected mean square mixes the
  random-slope and residual components with different coefficients), so
  the test aggregates – exactly as the two-level categorical engine
  reduces to per-subject difference vectors – to one within-group slope
  vector per group, then tests the equally weighted mean slope against
  the between-group slope variation by whole-group sign flips, with a
  parametric `F(1, ng - 1)` p-value for Euclidean dissimilarities. For a
  balanced design this is the exact two-stage slope test, validated to
  reproduce it to machine precision and to agree closely with
  `lmerTest`. The earlier attempt to build this denominator with
  [`structured_denominator()`](https://marine-ecologist.github.io/distmer/dev/reference/structured_denominator.md)
  was abandoned: the differing expected-mean-square coefficients
  inflated the synthesised denominator and pushed the pseudo-F ~1.6x
  below both gold standards, whereas the aggregated slope-vector test
  matches them exactly.
- [`permanova_random_slope()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_random_slope.md)
  (the balanced two-level categorical sign-flip engine) now rejects a
  continuous slope variable with a message pointing to the continuous
  engine, instead of failing obscurely.

*Internal checkpoint label: `2.6.0.9000` (unreleased).*

### Monte-Carlo asymptotic p-value: P(MC)

- [`permanova()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova.md)
  now reports `P(MC)` alongside `P(perm)` for designs tested in the full
  observation space (fixed, balanced mixed, nested, crossed, hybrid).
  Following McArdle & Anderson (2001) and Anderson & Robinson (2003),
  the numerator and denominator sums of squares are drawn from the
  eigenvalue- weighted chi-square distribution implied by the Gower
  matrix, so the p-value stays smooth when the true permutation space is
  tiny – the usual situation for the upper strata of a mixed design,
  where restricted permutation yields very few unique arrangements. It
  reduces exactly to the classical F distribution in the univariate
  Euclidean limit; validated against
  [`aov()`](https://rdrr.io/r/stats/aov.html) for a one-way design and
  against the exact balanced-mixed F-test over its stratum. Controlled
  by `mc_permutations` (default 9999; set 0 to disable) and reproducible
  from `seed`. Aggregated-unit engines (random slope, unequal random
  intercept, unbalanced nested/hybrid) test in a reduced centroid space
  where the full-matrix asymptotics do not hold, so their `P(MC)` is
  left `NA`.

### PRIMER-style nested containment labels

- Covariate-adjusted nested summaries now label strata with PRIMER
  containment notation (`Site(Location)`, `Area(Site(Location))`,
  `Residual`) in the table, the denominators, and the `V()` variance
  components, derived from the design’s nesting rather than hard-coded.
  The estimates of components of variation are reported in a single
  PRIMER-style block with `S()` fixed-effect sizes and `V()` random
  variance components.

*Internal checkpoint label: `2.5.0.9000` (unreleased).*

### Structured mixed-model denominators

- Added
  [`structured_denominator()`](https://marine-ecologist.github.io/distmer/dev/reference/structured_denominator.md),
  a Cornfield–Tukey / Satterthwaite denominator constructor that
  synthesises mean squares in containment-lattice order rather than by
  the minimum-norm match used in
  [`quasi_denominator()`](https://marine-ecologist.github.io/distmer/dev/reference/quasi_denominator.md).
  On a clean balanced design it returns the exact single-source
  denominator; on a perturbed (e.g. random-slope) design it returns the
  structured combination with Satterthwaite degrees of freedom.
  Validated against the clean balanced limit and against `lmerTest`
  denominator df (the latter requires the covariate to be
  within-centred, e.g. via
  [`decompose_covariate()`](https://marine-ecologist.github.io/distmer/dev/reference/decompose_covariate.md),
  so the fixed slope’s expected mean square is not contaminated by the
  random intercept).

### Covariate-adjusted nested inference

- Added
  [`permanova_nested_covariate()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_nested_covariate.md)
  (and auto-routing from
  [`permanova()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova.md)):
  categorical factors in a balanced `/`-generated hierarchy can now be
  tested while adjusting for one or more continuous covariates, matching
  PRIMER’s covariate handling (test the factors “above” the covariate
  relationship). Each factor’s numerator is adjusted for the covariate;
  its denominator stays the exact nested expected-mean-square stratum.
  The covariate is an observation-level nuisance charged to the deepest
  within-unit residual and tested there, so degrees of freedom partition
  every observation. A covariate with strong between-unit structure
  should be split into explicit between- and within-unit predictors. The
  [`permanova_nested()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_nested.md)
  error message now names a within-unit continuous covariate as the
  cause instead of advising `contr.sum`.
- Covariate-by-factor interactions are now tested. Writing
  `factor * covariate` in a nested design adds the homogeneity-of-slopes
  term (e.g. `Location:Volume`), tested over the residual as an
  observation-level term, with degrees of freedom taken from the
  residual. Previously such a term was silently dropped and the model
  collapsed to the additive one. Covariate-by- *random* interactions
  (random slopes) remain unsupported: their expected mean squares admit
  no exact or cleanly constructed denominator.
- Added
  [`decompose_covariate()`](https://marine-ecologist.github.io/distmer/dev/reference/decompose_covariate.md),
  a Mundlak-style helper that splits a covariate into between-unit
  (grouping-unit mean) and within-unit (deviation) parts, so a covariate
  with strong between-unit structure can be entered as two explicit
  predictors. The within part is a genuine observation-level covariate
  for
  [`permanova_nested_covariate()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_nested_covariate.md);
  the between part is a unit-level summary.

### Unified front door

- [`permanova()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova.md)
  accepts the dissimilarity on the formula’s left-hand side,
  `permanova(bray ~ treatment + (1 | site), data)`, matching the
  [`vegan::adonis2()`](https://vegandevs.github.io/vegan/reference/adonis.html)
  idiom; the one-sided `~ terms` form with a separate `distance`
  argument still works. The assembled PERMANOVA table is exposed as
  `$table` on the result.
- Added
  [`permanova()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova.md),
  a single entry point that accepts a formula (with data and distance)
  or a compiled `permanova_design`, routes it to the engine matching its
  random-effect provenance and balance, and returns a result that prints
  as a PERMANOVA table with a Factors block, denominators, and
  components of variation. Designs outside every supported envelope stop
  with a message rather than being approximated.
- Added `type = "I"|"II"|"III"` to select the fixed-effect
  sum-of-squares partition; it is reported as not applicable for designs
  whose denominator is fixed by the random structure. When called with a
  formula,
  [`permanova()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova.md)
  applies sum-to-zero contrasts to the fixed factors on a copy of the
  data, so Type III (and cross-level factor interactions) give the
  correct marginal hypotheses without the caller setting
  [`contrasts()`](https://rdrr.io/r/stats/contrasts.html). A
  pre-compiled design keeps whatever contrasts it was built with.
- The PERMANOVA table shows a grand `Total` row only when the reported
  sources partition every observation degree of freedom; aggregated-unit
  engines (random slope, unequal random intercept, unbalanced
  nested/hybrid) omit it.
- [`summary()`](https://rdrr.io/r/base/summary.html) of a fixed
  crossed/factorial model now prints the PRIMER-style
  expected-mean-squares, pseudo-F construction, and
  components-of-variation blocks. Each term’s EMS coefficient is the
  trace-derived `tr(H Z Z') / df` and its component estimate is
  `(MS_term - MS_res) / coefficient`; negative estimates are retained.
- Type II sums of squares now respect marginality instead of rejecting a
  main effect contained in an interaction. A term is adjusted for every
  term that does not contain it (the containing interactions are
  excluded), giving the standard SAS/`car`/PRIMER “Type II
  (conditional)” hypothesis; factor SS match `car::Anova(type = 2)` to
  numerical tolerance. No additive refit is required.
- Consolidated the vignette set from 24 example-indexed articles to 8
  organised by design family.

*Internal checkpoint label: `2.4.0.9000` (unreleased).*

### Development

- Preserved the sealed v2.0.0 technical candidate while preparing the
  publication handoff and the next validated inference phase.
- Closed Phase 23 external adversarial validation across Phases 19–25.
  The sealed v2.0.0 and Phase 25 artifacts matched their pinned
  checksums, all independent scalar/permutation references agreed to
  machine precision, and the confounded-strata safeguards held.
- Added an SVD condition-ratio guard for covariate-adjusted contrasts. A
  numerically near-aliased fixed model now rejects with an actionable
  package diagnostic before the contrast projector can leak a low-level
  [`solve()`](https://rdrr.io/r/base/solve.html) error; the ratio is
  retained in adjusted-mean metadata.
- Phase 25 targets child- and observation-level treatment inference in
  one unbalanced parent/child hierarchy with explicit unit-scale
  denominators and within-parent exchangeability contracts.
- Added
  [`permanova_unbalanced_nested_lower()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_unbalanced_nested_lower.md)
  with automatic assignment-level detection, connected-incidence and
  absolute-scale span guards, exhaustive or Monte Carlo restricted
  residual permutations, and complete Test Specifications.
- Extended mixed planned/pairwise contrasts and mixed global-null max-F
  to the same child- and observation-level unit models.
- Added scalar GLH and direct permutation reconstruction, a pinned
  two-level simulation audit, adversarial failure cases, and a Pyrifos
  occasion vignette.
- Phase 26 adds unbalanced treatment-by-region/site hybrid inference on
  equal-site-weight complete profile centroids, with additive region
  adjustment and region-by-exact-topology profile permutations.
- Routed mixed planned/pairwise and weak complete-null max-F inference
  through the same unbalanced profile model and common region reference
  grid.
- Added connected-incidence, complete-site-profile, topology-span,
  residual-df, zero-denominator, and non-Euclidean safeguards plus
  scalar and exhaustive profile reconstruction tests.
- Phase 27 adds the narrow correlated `(1 + condition | subject)` model
  for balanced complete two-condition repeated measures.
- Added residual-corrected random intercept–slope trace covariance
  estimates, optional positive-semidefinite truncation, and exhaustive
  whole-subject sign-flip inference under an explicit central-symmetry
  assumption.
- Added paired-t-squared, scalar covariance, direct sign-flip,
  asymmetric-null sensitivity, semimetric, and Pyrifos early-versus-late
  validation.
- Phase 28 adds
  [`permanova_pairwise_strong()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_pairwise_strong.md)
  for fixed one-factor all-pair complete-distribution equality
  hypotheses.
- Every intersection uses a separate sample-size-preserving raw-label
  action within connected equality components; closed testing then
  provides strong FWER without assuming subset pivotality.
- The proof excludes centroid-only nulls with unequal distributions,
  planned non-pairwise contrasts, dependent observations, and mixed or
  repeated designs. Exact reconstruction, partial-null preservation,
  calibration, and a Dune management vignette enforce that boundary.

*Internal checkpoint label: `2.0.0` (unreleased).*

### Validated support-contract consolidation

- Added
  [`permanova_capabilities()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_capabilities.md)
  and `SUPPORT-MATRIX.md` as machine- and human-readable maps of every
  primary inferential engine.
- Each registry row states its validated design envelope, estimand,
  denominator, permutation contract, error-control claim, and explicit
  exclusions; runtime design validators remain authoritative.
- Extended deterministic and structural release gates to cover the
  registry, support matrix, public API, all 20 vignettes, and the
  complete Phase 1–22 documentation sequence.

*Internal checkpoint label: `1.8.0` (unreleased).*

### Mixed-design global-null max-F inference

- Added
  [`permanova_mixed_maxT()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_mixed_maxT.md)
  and
  [`permanova_mixed_maxT_pairwise()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_mixed_maxT_pairwise.md)
  for all balanced and unbalanced mixed treatment structures already
  supported by the mixed-contrast engine.
- One shared raw experimental-unit permutation stream respects
  whole-unit, exact-size, hierarchy-topology, block-by-size, or region
  strata and forms a maximum pseudo-F distribution.
- Results explicitly claim only weak familywise error control under the
  complete treatment-null randomization hypothesis.
- Added exact block enumeration, scalar unit-centroid reconstruction,
  non-Euclidean coverage, a 5,000-dataset weak-FWER audit, and a Pyrifos
  dose vignette.

*Internal checkpoint label: `1.7.0` (unreleased).*

### Explicit continuous-covariate adjusted contrasts

- Added
  [`derive_adjusted_means()`](https://marine-ecologist.github.io/distmer/dev/reference/derive_adjusted_means.md),
  [`permanova_adjusted_contrasts()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_adjusted_contrasts.md),
  and
  [`permanova_adjusted_pairwise()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_adjusted_pairwise.md)
  for full-rank fixed models containing categorical factors and ordinary
  numeric covariates.
- Every continuous nuisance covariate requires an explicit finite
  reference value; categorical nuisance factors use declared equal-grid
  or observed- frequency weighting.
- Added coefficient-space GLH checks, exhaustive direct Freedman–Lane
  reconstruction, non-Euclidean coverage, a pinned 500-replicate
  validation, and a Dune Management-by-A1 vignette.

*Internal checkpoint label: `1.6.0` (unreleased).*

### Unbalanced mixed-design contrasts

- Extended
  [`permanova_mixed_contrasts()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_mixed_contrasts.md)
  and
  [`permanova_mixed_pairwise()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_mixed_pairwise.md)
  to unequal random-intercept units, unbalanced top-level parent/child
  hierarchies, and incomplete or unequally replicated randomized blocks.
- Cluster-size, complete-topology, and block-by-cell-size strata reuse
  the amended absolute-scale span guard before any contrast permutation.
- Incomplete blocks use equally weighted observed cell centroids, an
  additive block model, and treatment means averaged over one common
  block reference grid.
- Added scalar GLH checks for all three structures, exhaustive direct
  permutation reconstruction, a pinned 500-replicate audit, and an
  unequal- replication Pyrifos vignette.

*Internal checkpoint label: `1.5.1` (unreleased).*

### Stratified-exchangeability safeguard amendment

- Replaced relative QR rank checks in the Phase 11 and Phase 12
  exchangeability-span guards with an SVD rank threshold measured
  absolutely relative to the original treatment-design scale.
- Fully size-confounded, hierarchy-topology-confounded, and within-block
  cell-size-confounded treatments now reject before permutation instead
  of admitting a degenerate p-value of one.
- Added adversarial floating-dust regression tests for all three
  validators; supported cluster-centroid, top-centroid, and
  cell-centroid references remain unchanged to machine precision.

*Internal checkpoint label: `1.5.0` (unreleased).*

### Release-candidate hardening

- Added multi-platform R CMD check CI for release R on Linux, macOS, and
  Windows, including dependency caching and source-package artifact
  upload.
- Added a structural release audit for version/news parity,
  public-API/Rd parity, vignette count, S3 registration, build hygiene,
  and placeholder publication metadata.
- Added release-readiness documentation separating a green technical
  package gate from the unresolved maintainer identity, repository URL,
  and bug-tracker metadata required before public submission.
- Added deterministic release-audit tests and a bounded benchmark smoke
  gate.

*Internal checkpoint label: `1.4.0` (unreleased).*

### Global-null max-F inference

- Added
  [`permanova_maxT()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_maxT.md)
  and
  [`permanova_maxT_pairwise()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_maxT_pairwise.md)
  for fixed one-factor planned or pairwise contrast families.
- Every contrast uses one shared raw-observation permutation stream;
  adjusted p-values compare the observed pseudo-F with the permutation
  maximum.
- Results explicitly label the guarantee as weak FWER under the complete
  equal-centroid null and do not claim strong subset-null control.
- Added exact brute-force joint-distribution checks, non-Euclidean
  coverage, a 5,000-dataset family-error audit, and a Dune Management
  vignette.

*Internal checkpoint label: `1.3.0` (unreleased).*

### Mixed-design planned contrasts

- Added
  [`permanova_mixed_contrasts()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_mixed_contrasts.md)
  and
  [`permanova_mixed_pairwise()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_mixed_pairwise.md)
  for one categorical treatment in four validated balanced structures: a
  random intercept, a top-level nested hierarchy, a declared randomized
  block, and a crossed-nested region/site hybrid.
- Contrasts operate on equally weighted cluster, top-unit,
  treatment-cell, or treatment-profile centroids and use the
  corresponding experimental-unit, treatment-by-block, or
  treatment-by-region residual denominator.
- Contrast-specific Freedman–Lane residual permutations move whole units
  or remain within blocks/regions. Unsupported imbalance and wrong-level
  treatments reject explicitly.
- Added scalar GLH references for every structure, exhaustive
  response-level permutation agreement, a pinned null/power validation,
  and a Pyrifos planned-contrast vignette.

*Internal checkpoint label: `1.2.0` (unreleased).*

### Multifactor marginal contrasts

- Added
  [`derive_marginal_means()`](https://marine-ecologist.github.io/distmer/dev/reference/derive_marginal_means.md)
  with explicit equal-reference-grid or observed-frequency nuisance
  weighting for full-rank categorical fixed models.
- Added
  [`permanova_marginal_contrasts()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_marginal_contrasts.md)
  and
  [`permanova_marginal_pairwise()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_marginal_pairwise.md)
  with general-linear-hypothesis projectors, contrast-specific
  Freedman–Lane permutations, and explicit multiplicity adjustment.
- Continuous nuisance variables, random effects, aliased reference
  grids, and implicit weighting reject before inference.
- Added scalar GLH and response-reconstruction references, a
  500-replicate null/power validation, and a Dune Management-by-Use
  vignette.

*Internal checkpoint label: `1.1.0` (unreleased).*

### Unequal-random-intercept calibration audit

- Added
  [`audit_unequal_random_intercept_calibration()`](https://marine-ecologist.github.io/distmer/dev/reference/audit_unequal_random_intercept_calibration.md)
  to compare exhaustive unique treatment-label assignments, Monte Carlo
  inference, and the complete distance engine across non-overlapping
  seed batches.
- Proved and regression-tested that the intercept-only cluster-centroid
  Freedman–Lane transformation is algebraically identical to raw
  centroid permutation.
- Added Gaussian, heavy-tail, unequal-allocation, and deliberate
  dispersion- violation scenarios. Only exchangeable nulls contribute to
  the calibration gate.
- Distinguished all profile rearrangements from unique treatment-label
  assignments in Phase 11 result metadata.
- Fixed a zero-permutation matrix-construction warning in stratified
  plans.
- The 40,000-dataset audit explains the earlier elevated 500-dataset
  batch: valid exact rejection rates are 4.53–4.76%, with full-engine
  agreement to machine precision.

*Internal checkpoint label: `1.0.0` (unreleased; not the future stable
release).*

### Planned contrasts and pairwise post-hoc inference

- Added
  [`permanova_contrasts()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_contrasts.md)
  for user-specified one-df contrasts among the centroids of one
  categorical fixed factor with unequal group sizes.
- Each contrast retains all orthogonal group-centroid directions as
  nuisance structure and receives its own Freedman–Lane residual
  permutation.
- Added `permanova_pairwise()` to construct every pairwise centroid
  difference. Holm is the default adjustment; Bonferroni, BH, BY, and no
  adjustment are explicit alternatives.
- Invalid/non-zero-sum contrast definitions, mixed or multifactor
  designs, aliases, and zero residual df reject before inference.
- Added exact constrained-scalar references, non-Euclidean coverage, a
  pinned null/power validation, and a Dune planned/pairwise vignette.

*Internal checkpoint label: `0.9.0` (unreleased).*

### Unbalanced nested and randomized-block inference

- Added
  [`permanova_unbalanced_nested()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_unbalanced_nested.md)
  for one treatment assigned to top units of an unequal parent/child
  hierarchy. Equally weighted top centroids exchange only among
  identical sorted child-count topologies.
- Added
  [`permanova_unbalanced_block()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_unbalanced_block.md)
  for incomplete and/or unequally replicated randomized blocks with an
  explicit treatment:block component. Equally weighted intact cell
  centroids exchange only within block and exact cell-size strata.
- Both engines require their restricted permutation action to span every
  treatment contrast and reject disconnected incidence, confounded
  topology, aliases, or zero residual df.
- Added exhaustive scalar reconstructions, non-Euclidean tests, a
  combined operating-characteristic harness, and a Dune unbalanced-mixed
  vignette.

*Internal checkpoint label: `0.8.0` (unreleased).*

### Unequal-replication random-intercept inference

- Added
  [`permanova_unequal_random_intercept()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_unequal_random_intercept.md)
  for one cluster-level categorical treatment and one random intercept
  with unequal observations per unit.
- The estimand is explicitly cluster weighted: the analysis uses equally
  weighted cluster centroids rather than silently giving larger clusters
  more influence.
- Freedman–Lane residual permutations move complete cluster centroids
  only within exact cluster-size strata; those strata must span every
  treatment contrast.
- No pooled observation-weighted EMS or parametric F distribution is
  claimed. Non-Euclidean empirical inference retains the complete Gower
  spectrum.
- Added exhaustive scalar references, unsafe-design failures, a Pyrifos
  vignette, and pinned null/power validation.

*Internal checkpoint label: `0.7.0` (unreleased).*

### Unbalanced fixed-effect inference

- Added
  [`derive_unbalanced_hypotheses()`](https://marine-ecologist.github.io/distmer/dev/reference/derive_unbalanced_hypotheses.md)
  and
  [`permanova_unbalanced()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_unbalanced.md)
  for explicit fully estimable Type I, II, or III fixed-effect
  hypotheses.
- Type I tests are sequential; additive Type II tests are marginal and
  order-invariant; Type III tests retain all other terms and require a
  full-rank model plus explicit sum-to-zero contrasts when interactions
  occur.
- Partially estimable terms, aliased incomplete Type III models, Type II
  main effects contained in interactions, zero residual df, and all
  random-effect structures reject before inference.
- Every term uses the common full-model residual denominator and its own
  Freedman–Lane reduced-model residual permutation.
- Added exact `5!` enumeration, non-Euclidean empirical inference,
  scalar Type I/II/III references, incomplete-additive coverage, and
  pinned unequal- cell null/power validation.

*Internal checkpoint label: `0.6.0` (unreleased).*

### Balanced crossed-nested hybrids

- Added
  [`permanova_hybrid()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_hybrid.md)
  for one fixed treatment crossed with regions and with sites nested
  within regions.
- The exact envelope requires explicit region, region:site,
  treatment:region, and treatment:region:site random-intercept
  covariance components, complete treatment occurrence within every
  site, and equal residual replication.
- Added orthogonal crossed-nested strata, trace-derived EMS, the exact
  `MS[treatment:region]` denominator, and unrestricted five-component
  variance estimates.
- Exact or seeded Monte Carlo permutations move complete
  treatment-by-region profiles within regions, carrying every site and
  replicate descendant intact.
- Added exhaustive profile-permutation references, non-Euclidean
  inference, strict provenance/balance failures, and pinned null/power
  validation.

*Internal checkpoint label: `0.5.0` (unreleased).*

### Design-aware multivariate dispersion

- Added
  [`permanova_dispersion()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_dispersion.md)
  for one categorical fixed target in balanced one-way,
  single-random-intercept, nested, and canonical randomized-block
  designs.
- Distances to group spatial medians are the default; centroids and the
  equal-group-size bias adjustment are explicit options.
- Positive and negative principal-coordinate axes are retained. Negative
  squared pseudo-Euclidean distances to a center are counted, warned
  about, and changed to zero rather than being silently discarded.
- Dispersion responses use the target term’s exact EMS denominator and
  full-model residual permutations at the correct exchangeability level.
- One-way geometry, F statistics, and permutation vectors agree with
  [`vegan::betadisper()`](https://vegandevs.github.io/vegan/reference/betadisper.html)
  and
  [`permutest()`](https://vegandevs.github.io/vegan/reference/anova.cca.html)
  under identical arrangements.
- Added exhaustive scalar-reference checks for one-way and crossed
  residual permutations plus a four-structure operating-characteristic
  validation.
- Improved the Phase 6 cross-level interaction diagnostic to name
  explicit `contr.sum` assignment and added an executable split-plot
  example.

*Internal checkpoint label: `0.4.0` (unreleased).*

### Crossed randomized-block permutation inference

- Added `p_method = "permutation"` to
  [`permanova_crossed()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_crossed.md)
  for the canonical balanced fixed-treatment × random-block model with
  an explicitly declared treatment-by-block random component.
- Permutations rearrange complete replicated treatment cells within each
  block; no replicate separates from its cell and no cell crosses a
  block.
- Added exact `(treatments!)^blocks` enumeration, seeded Monte Carlo
  sampling, attainable-space metadata, duplicate-draw reporting, and
  retained statistics.
- Empirical p-values remain available for non-Euclidean dissimilarities
  while the complete negative spectrum is retained.
- Added exhaustive two- and three-treatment references, direct
  response-level statistic comparisons, exchangeability failures, a
  pinned validation with dispersion sensitivity, and a randomized-block
  ecological vignette.

*Internal checkpoint label: `0.3.0` (unreleased).*

### Level-specific balanced nested inference

- Generalized
  [`permanova_nested()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_nested.md)
  from top-level effects to fixed terms assigned at any one orthogonal
  level of a balanced hierarchy.
- Region effects permute whole regions, site effects permute intact
  sites within regions, and observation effects permute observations
  within deepest units.
- Test Specifications now record assignment level, parent level, exact
  level-specific EMS denominator, attainable permutation space, and the
  hypothesis-specific exchangeability algorithm.
- Added exhaustive global, within-parent, and within-unit permutation
  tests; classical nested-ANOVA comparisons at all three levels;
  balanced cross-level interaction checks; and explicit rejection of
  level-spanning covariates.
- Added
  [`simulate_nested_level_validation()`](https://marine-ecologist.github.io/distmer/dev/reference/simulate_nested_level_validation.md)
  with binomial uncertainty and per-level classical agreement
  diagnostics.

*Internal checkpoint label: `0.2.0` (unreleased).*

### Balanced nested mixed designs

- Added exact inference for one complete balanced `/`-generated
  random-intercept hierarchy, such as `(1 | region/site)`, when fixed
  predictors are constant within top-level units.
- Added trace-derived nested EMS, an exact residual top-unit
  denominator, and term-specific Freedman–Lane permutation of intact
  top-level hierarchies.
- Added unrestricted nested variance-component estimates, with zero
  truncation available only as a separately reported sensitivity
  analysis.
- Unsupported imbalance, incomplete or crossed hierarchies, random
  slopes, and wrong-level fixed predictors reject explicitly.
- Added analytical ANOVA agreement, exhaustive `6!` permutation
  coverage, three-level hierarchy checks, a 500-replicate calibration
  experiment, and a nested ecological-community vignette.

*Internal checkpoint label: `0.1.1` (unreleased).*

### Crossed-EMS safety amendment

- Explicit `:` random interactions are distinct from `/` nesting
  throughout the canonical Design Object.
- Declared fixed-by-random interaction components enter trace-derived
  EMS and can provide exact classical interaction denominators.
- Omitted fixed-by-random components now reject by default. Additive
  pooling requires `assume_additive = TRUE` and records all assumed-zero
  components and pooled strata.
- Satterthwaite p-values are suppressed for non-Euclidean
  dissimilarities.
- Added deterministic adversarial checks and a reproducible
  crossed-interaction Type I validation with binomial confidence
  intervals.
