# distmer

A transparent, extensible R framework for distance-based fixed and mixed models
(PERMANOVA-style) in ecological community analysis. Unsupported inference
structures fail explicitly instead of receiving an implicit denominator or
permutation scheme.

## Quick start

`permanova()` is the main user interface. Give it an `lmer`-style formula with
the dissimilarity on the left-hand side. It validates the sampling structure,
matches the model to a supported distance-based inference method, and returns a
familiar results table with a permutation p-value and, where applicable, a
smooth Monte-Carlo `P(MC)`.

```r
library(distmer)
data(dune, dune.env, package = "vegan")
bray <- vegan::vegdist(dune, method = "bray")

# one fixed factor (agrees with vegan::adonis2)
permanova(bray ~ Management, dune.env, permutations = 999, seed = 1)

# unbalanced fixed factors, Type II sums of squares
permanova(bray ~ Management + Moisture, dune.env, type = "II")

# one random intercept (whole-plot EMS denominator)
permanova(bray ~ treatment + (1 | site), sampling_design)

# split-plot / RCB ANCOVA (covariate adjusted)
permanova(bray ~ treatment + volume + (1 | site), sampling_design)

# nested hierarchy, covariate-adjusted nested strata, crossed/randomized block,
# crossed-nested hybrid, random slopes -- all through the same call:
permanova(bray ~ treatment + (1 | region/site), sampling_design)
permanova(bray ~ x + (1 + x | group), sampling_design)
```

The selected implementation is printed in the header and stored in
`fit$engine`; its detailed result is in `fit$fit`, and `summary(fit)` adds
estimates of components of variation. Designs outside every supported envelope
stop with a message naming the reason rather than guessing a denominator.

Advanced users can compile a design explicitly with
`compile_permanova_design()` and call a specialised function when they need a
low-level argument. For example, `permanova_dispersion()` and the contrast
families take a compiled design. In day-to-day use, start with `permanova()`.

## Balanced random-intercept inference

```r
fit <- permanova(bray ~ treatment + (1 | site), sampling_design,
                 permutations = 9999, seed = 5)
fit$fit$tests$treatment$symbolic_ems
```

Use `permutations = "all"` for exhaustive enumeration when the whole-plot
permutation space does not exceed `max_exact`.

## Balanced crossed mixed inference

Phase 3 supports complete, equally replicated factorial designs whose
fixed-by-random covariance components are declared explicitly:

```r
design <- compile_permanova_design(
  ~ treatment +
    (1 | block) + (1 | observer) +
    (1 | treatment:block) + (1 | treatment:observer),
  data = sampling_design,
  distance = community_distance
)

fit <- permanova_crossed(design)
fit$variance_components$unrestricted
```

Trace-derived EMS coefficients select and pool exact orthogonal denominator
strata. `quasi_denominator()` constructs an EMS-matched linear combination when
no single exact stratum exists. Negative unrestricted variance components are
retained; request `truncate = TRUE` to report a separate zero-truncated
sensitivity result. Parametric Satterthwaite p-values are never implicit and
require `p_method = "satterthwaite"`; they are suppressed for non-Euclidean
dissimilarities. Omitted fixed-by-random components fail by default. The
`assume_additive = TRUE` escape hatch is an explicit zero-variance assumption,
recorded with the pooled denominator strata in every Test Specification.

The canonical randomized-block model also has a permutation path:

```r
design <- compile_permanova_design(
  ~ treatment + (1 | block) + (1 | treatment:block),
  data = sampling_design,
  distance = community_distance
)

fit <- permanova_crossed(
  design,
  p_method = "permutation",
  permutations = 9999,
  seed = 5
)
```

Complete replicated treatment cells move among treatment positions within each
block. Use `permutations = "all"` when `(treatments!)^blocks` does not exceed
`max_exact`. This empirical path supports retained non-Euclidean distances.
Additional fixed/random factors, nesting, imbalance, or additive omissions
remain explicit failures for permutation inference.

## Balanced nested mixed inference

```r
fit <- permanova(
  community_distance ~ region_treatment + site_treatment +
    observation_treatment + (1 | region/site),
  sampling_design, permutations = 9999, seed = 5
)
fit$fit$assignment_levels
fit$fit$tests$site_treatment$symbolic_ems
fit$fit$variance_components$unrestricted
```

`derive_nested_ems()` exposes the trace-derived numerator and denominator EMS.
`estimate_nested_variance_components()` solves the hierarchy's balanced EMS
system and retains negative unrestricted estimates; `truncate = TRUE` adds a
separate zero-truncated sensitivity result.

Every term's estimable projector must lie wholly in one stratum. Balanced
cross-level interactions are supported when parameterized with orthogonal
contrasts such as explicitly assigned `contr.sum`; otherwise the engine rejects
the ambiguous projector and names that remedy. Contrasts are never changed
silently because doing so could change the tested hypothesis.

## Balanced crossed-nested hybrid inference

Phase 9 combines one nested hierarchy with explicitly declared treatment
interactions:

```r
design <- compile_permanova_design(
  ~ treatment +
    (1 | region/site) +
    (1 | treatment:region) +
    (1 | treatment:region:site),
  data = sampling_design,
  distance = community_distance
)

hybrid <- permanova_hybrid(
  design,
  permutations = 9999,
  seed = 5
)
```

The treatment numerator is tested over `MS[treatment:region]`. Permutations
move complete treatment-by-region profiles within each region, carrying every
nested site and replicate intact. The first hybrid envelope is deliberately
exact: one treatment, one region/site hierarchy, complete treatment occurrence
within every site, equal replication, and all treatment-by-hierarchy covariance
components declared explicitly.

## Unbalanced fixed-effect inference

Phase 10 makes the requested hypothesis explicit:

```r
fit <- permanova(
  community_distance ~ habitat + season,
  sampling_design,
  type = "II",
  permutations = 9999,
  seed = 5
)
```

Type I tests follow formula order. Type II tests respect marginality: a main
effect is adjusted for terms that do not contain it while its containing
higher-order interactions are held out of that reduced model. Type III retains
all other terms and requires a full-rank model; with interactions, assign
explicit sum-to-zero contrasts such as `contr.sum` before compiling.
`derive_unbalanced_hypotheses()` exposes df, nominal df, adjustment terms,
aliases, and the projector before inference.

Every selected term must be fully estimable and is tested over the final
full-model residual mean square with its own Freedman--Lane permutation. Phase
10 does not infer an unbalanced mixed-model EMS denominator: any random effect
causes an explicit failure.

## Unequal-replication random-intercept inference

Phase 11 targets equally weighted whole-plot centroids when observation counts
differ among random-intercept units:

```r
design <- compile_permanova_design(
  ~ treatment + (1 | site),
  data = sampling_design,
  distance = community_distance
)

fit <- permanova_unequal_random_intercept(
  design,
  permutations = 9999,
  seed = 5
)
```

The pseudo-F uses a cluster-centroid residual denominator. Freedman--Lane
residuals move only among whole clusters with the same observation count, and
the size-stratified permutation action must span every treatment contrast.
The result is a cluster-weighted empirical test. It is not an
observation-weighted EMS approximation, and no parametric F p-value is
reported.

## Unbalanced nested and randomized-block inference

Phase 12 uses separate empirical engines for two different unit structures:

```r
nested_fit <- permanova_unbalanced_nested(
  compile_permanova_design(
    ~ treatment + (1 | region/site), sampling_design, community_distance
  ),
  permutations = 9999,
  seed = 5
)

block_fit <- permanova_unbalanced_block(
  compile_permanova_design(
    ~ treatment + (1 | block) + (1 | treatment:block),
    sampling_design, community_distance
  ),
  permutations = 9999,
  seed = 5
)
```

The nested engine gives equal weight to top-unit centroids and permutes
residuals only among identical sorted child observation-count topologies. The
block engine gives equal weight to observed treatment-by-block cell centroids
and permutes intact profiles only within block and exact cell size. In both
cases the restricted permutation action must span every treatment contrast;
otherwise the global hypothesis rejects before inference.

## Lower-level unbalanced nested inference

Phase 25 handles a categorical treatment assigned below the parent in one
unbalanced parent/child hierarchy:

```r
lower_fit <- permanova_unbalanced_nested_lower(
  compile_permanova_design(
    ~ treatment + (1 | region/site), sampling_design, community_distance
  ),
  level = "auto",
  permutations = 9999,
  seed = 5
)
```

If treatment is constant within child, equally weighted child centroids are
adjusted for parent and residuals move only among children in the same parent
and exact-size stratum. If treatment varies within child, observations are
adjusted for child and residuals move only within child. The same unit models
are available to `permanova_mixed_contrasts()`,
`permanova_mixed_pairwise()`, and the explicitly weak
`permanova_mixed_maxT()` family. Parent-level treatment remains the existing
topology-stratified top-unit problem.

## Unbalanced crossed-nested hybrid inference

Phase 26 handles one treatment crossed with regions and sites nested within
regions when site counts, profile incidence, or cell replication are unequal:

```r
hybrid_fit <- permanova_unbalanced_hybrid(
  compile_permanova_design(
    ~ treatment + (1 | region/site) +
      (1 | treatment:region) + (1 | treatment:region:site),
    sampling_design, community_distance
  ),
  permutations = 9999,
  seed = 5
)
```

Each observed treatment-by-region profile must contain every site in that
region. Observations are first averaged within treatment-by-site cells, then
site-cell centroids receive equal weight in the profile. The treatment test is
adjusted for region and uses the remaining profile-centroid mean square;
inference is empirical and makes no unbalanced EMS claim. Intact profiles move
only within region and the same ordered site-specific cell-count topology.
Missing whole profiles are allowed when incidence remains connected and those
restricted actions span every treatment contrast. Mixed planned/pairwise and
weak complete-null max-F contrasts reuse the same profiles and a common region
reference grid.

## Correlated random intercepts and slopes

Phase 27 supports a balanced complete two-condition repeated-measures model:

```r
slope_fit <- permanova(
  community_distance ~ condition + (1 + condition | subject),
  repeated_design,
  permutations = "all",
  truncate = TRUE
)
```

A continuous covariate slope, `~ x + (1 + x | group)`, routes instead to
`permanova_random_slope_continuous()`, which aggregates to per-group slope
vectors and reproduces the exact two-stage slope test.

The random intercept is the subject midpoint under centred two-condition
coding; the random slope is the complete second-minus-first condition
difference vector. The fixed-effect pseudo-F tests the population mean
difference over between-subject slope variation and equals a paired t statistic
squared for a Euclidean scalar response. Inference sign-flips whole subject
difference vectors and therefore requires central symmetry under the null; it
is not described as unrestricted condition-label exchangeability. The reported
trace covariance removes equal within-cell sampling noise by method of moments,
retains negative unrestricted estimates, and optionally supplies a separate
positive-semidefinite sensitivity.

More than two conditions, incomplete subjects, unequal or unreplicated cells,
extra fixed effects, and asymmetric difference distributions lie outside this
contract. Non-Euclidean empirical inference is retained, but its Gaussian
parametric p-value is suppressed.

## Planned contrasts and pairwise comparisons

Phase 13 tests explicit one-df contrasts among the centroids of one categorical
fixed factor, including unequal group sizes:

```r
planned <- rbind(
  `A - B` = c(A = 1, B = -1, C = 0),
  `AB - C` = c(A = 0.5, B = 0.5, C = -1)
)

planned_fit <- permanova_contrasts(
  design,
  contrasts = planned,
  adjust = "holm",
  permutations = 9999,
  seed = 5
)

pairwise_fit <- permanova_pairwise(
  design,
  adjust = "holm",
  permutations = 9999,
  seed = 5
)
```

Contrast coefficients apply to unweighted group centroids and must sum to
zero. Each contrast receives a reduced model that retains every orthogonal
group-centroid direction plus its own Freedman--Lane permutation. Holm is the
default multiplicity adjustment; Bonferroni, BH, BY, and no adjustment are
available explicitly.

## Multifactor marginal contrasts

Phase 15 defines how nuisance factors are averaged before testing a target-
factor contrast:

```r
means <- derive_marginal_means(
  design, factor = "treatment", weights = "equal"
)

fit <- permanova_marginal_contrasts(
  design,
  factor = "treatment",
  contrasts = planned,
  weights = "equal",
  permutations = 9999,
  seed = 5
)

pairs <- permanova_marginal_pairwise(
  design,
  factor = "treatment",
  weights = "proportional",
  permutations = 9999,
  seed = 5
)
```

Equal weights average the complete categorical nuisance grid. Proportional
weights use observed joint nuisance-cell frequencies. The fixed model must be
full rank and categorical. Random effects reject; use the explicit covariate-
reference interface below when numeric predictors are present.

## Continuous-covariate adjusted contrasts

Phase 20 tests fixed-model contrasts at user-declared numeric covariate values:

```r
means <- derive_adjusted_means(
  design, factor = "treatment", at = c(temperature = 20),
  weights = "equal"
)

fit <- permanova_adjusted_contrasts(
  design, factor = "treatment", contrasts = planned,
  at = c(temperature = 20), weights = "equal",
  permutations = 9999, seed = 5
)

pairs <- permanova_adjusted_pairwise(
  design, factor = "treatment", at = c(temperature = 20),
  weights = "proportional", permutations = 9999, seed = 5
)
```

Every ordinary numeric nuisance covariate requires one named finite scalar in
`at`; none is inferred from the observations. Categorical nuisance factors
are averaged over an equal reference grid or observed joint frequencies.
Random effects, aliased or ill-conditioned fixed models, unsupported predictor
classes, and incomplete reference specifications reject explicitly. The
adjusted-mean result reports `fixed_model_condition_ratio`; models whose ratio
cannot support stable information-matrix construction at `tol` reject before
projector construction.

## Mixed-design contrasts

Phase 16 tests treatment contrasts on the experimental-unit scale selected by
the declared balanced random structure:

```r
fit <- permanova_mixed_contrasts(
  design,
  contrasts = planned,
  adjust = "holm",
  permutations = 9999,
  seed = 5
)

pairs <- permanova_mixed_pairwise(
  design, adjust = "holm", permutations = 9999, seed = 5
)
```

One random intercept uses equally weighted whole-unit centroids; a top-level
nested treatment uses intact top-unit centroids; a declared randomized block
uses treatment-by-block cell centroids; and the balanced region/site hybrid
uses treatment-by-region profiles. Unequal unit sizes restrict whole-unit
residuals by exact size, unbalanced nested units by complete topology, and
incomplete block cells jointly by block and exact cell size. Incomplete blocks
use one common block reference grid for treatment means. Treatments assigned
below the top of an unbalanced hierarchy remain explicit failures.

## Global-null max-F contrasts

Phase 17 uses one shared raw-observation permutation distribution for a fixed
one-factor contrast family:

```r
fit <- permanova_maxT(
  design, contrasts = planned, permutations = 9999, seed = 5
)

pairs <- permanova_maxT_pairwise(
  design, permutations = 9999, seed = 5
)
```

Each permutation contributes the largest contrast pseudo-F. `p_maxT` compares
each observed statistic with those maxima. The guarantee is deliberately
limited to weak familywise error control under the complete equal-centroid
null. Strong subset-null control is not claimed by this engine.

## Strong closed pairwise inference

Phase 28 provides a narrower family for which subset-null exchangeability can
be proved:

```r
pairs <- permanova_pairwise_strong(
  design, permutations = 9999, seed = 5
)
```

For every intersection of pairwise equalities, factor levels joined by those
equalities form connected components. Observations are reassigned only within
each component while preserving group sizes; non-null groups remain fixed.
Under independent observations and equality of the complete distributions
inside each component, this raw-label action leaves the intersection-null law
invariant. A maximum pseudo-F tests the intersection, and closed testing gives
`p_closed` strong FWER without invoking subset pivotality.

The distinction is substantive: this is not a strong-control claim for equal
centroids when dispersions or distribution shapes differ. Planned non-pairwise
contrasts, dependent observations, repeated measures, and mixed designs remain
unsupported by the strong engine.

## Mixed-design global-null max-F contrasts

Phase 21 carries the same complete-null maximum-statistic contract to the
experimental-unit scale of supported mixed designs:

```r
fit <- permanova_mixed_maxT(
  design, contrasts = planned, permutations = 9999, seed = 5
)

pairs <- permanova_mixed_maxT_pairwise(
  design, permutations = 9999, seed = 5
)
```

One shared raw experimental-unit permutation stream is used across every
contrast. Balanced whole units move intact; unequal units remain within exact
size strata; unbalanced nested top units remain within topology strata; and
block or hybrid profiles remain within the corresponding block, block-by-size,
or region strata. `p_maxT` provides weak FWER under the complete treatment-null
randomization hypothesis. It does not provide strong subset-null control.

## V2 support registry

Phase 22 consolidates the public inference contract without broadening any
runtime support envelope:

```r
permanova_capabilities()
```

Each row names a primary engine and its validated design, estimand,
denominator, permutation contract, error-control claim, and explicit
exclusions. The human-readable companion is
`inst/SUPPORT-MATRIX.md`. These maps aid discovery; the compiler and each
engine's runtime validator remain authoritative for a particular design.

## Calibration hardening

Phase 14 audits unequal-random-intercept inference against the exact unique
treatment-label randomization distribution across independent seed batches:

```r
calibration_audit <- audit_unequal_random_intercept_calibration(
  batches = 20,
  simulations_per_batch = 500,
  permutations = 99,
  engine_checks = 10,
  seed = 20260719
)

marginal_validation <- simulate_marginal_contrast_validation(
  simulations = 500,
  permutations = 99,
  seed = 20260719
)
```

## Release readiness

The source tree is the 0.9.0.9000 pre-release development version. The next
stable version is 1.0.0; no tag or submission is implied by development checks.
See `inst/RELEASE-HANDOFF.md` for the owner-only publication boundary.
The package version policy is recorded in `inst/VERSIONING.md`.

The structural audit, test suite, source build, and package check are
reproducible locally:

```sh
Rscript tools/release-audit.R
Rscript -e 'testthat::test_local(reporter = "progress")'
R CMD build .
R CMD check --no-manual distmer_0.9.0.9000.tar.gz
```

Hosted CI repeats an `--as-cran` check on Linux, macOS, and Windows. The audit
reports publication metadata separately from technical package integrity; set
`STRICT_RELEASE_METADATA=true` to make unresolved metadata fatal.
See `inst/RELEASE-READINESS.md` for the exact handoff gate.

The current gate covers the deterministic validation battery, complete source
test suite, all seven consolidated vignettes, source build, and installed
package check. Results are recorded only after those commands have run against
the active version; retired Phase artifacts are not treated as evidence for
the development tree.

The audit separates `13,824` distinct profile rearrangements from `216` unique
labels in the canonical balanced-allocation example, verifies identical exact
p-values, and records seed ranges spaced one million apart. Exchangeable
Gaussian, heavy-tail, and unequal-allocation nulls determine the calibration
gate. Treatment-specific dispersion is reported only as a violated-
exchangeability sensitivity.

## Design-aware multivariate dispersion

`permanova_dispersion()` tests whether mean distances to group spatial medians
or centroids differ while retaining the design's exact EMS denominator and
restricted exchangeability units:

```r
design <- compile_permanova_design(
  ~ treatment + (1 | block) + (1 | treatment:block),
  data = sampling_design,
  distance = community_distance
)

dispersion <- permanova_dispersion(
  design,
  term = "treatment",
  type = "median",
  permutations = 9999,
  seed = 5
)
```

The original positive and negative principal-coordinate axes are retained.
The default spatial median has materially better small-sample null calibration
than the centroid in the pinned validation; centroid results remain available
as an explicit sensitivity. Phase 8 requires a design containing exactly one
categorical fixed main effect with equal target-group sizes. A significant
dispersion result can complicate interpretation of a location test, while a
non-significant result does not prove exchangeability.

## Validation, simulation, and benchmarks

```r
run_permanova_validation()

simulation <- simulate_mixed_validation(
  simulations = 1000,
  permutations = 999,
  seed = 20260718
)

interaction_validation <- simulate_crossed_interaction_validation(
  simulations = 500,
  blocks = 8,
  replicates = 3,
  seed = 20260718
)

nested_validation <- simulate_nested_validation(
  simulations = 500,
  regions_per_treatment = 4,
  sites_per_region = 2,
  replicates = 3,
  seed = 20260718
)

level_validation <- simulate_nested_level_validation(
  simulations = 500,
  regions_per_treatment = 4,
  replicates = 3,
  seed = 20260718
)

crossed_permutation_validation <- simulate_crossed_permutation_validation(
  simulations = 500,
  permutations = 99,
  blocks = 8,
  replicates = 3,
  seed = 20260718
)

unequal_random_intercept_validation <-
  simulate_unequal_random_intercept_validation(
    simulations = 500,
    permutations = 99,
    seed = 20260719
  )

unbalanced_mixed_validation <- simulate_unbalanced_mixed_validation(
  simulations = 500,
  permutations = 99,
  seed = 20260719
)

contrast_validation <- simulate_contrast_validation(
  simulations = 500,
  permutations = 99,
  seed = 20260719
)

calibration_audit <- audit_unequal_random_intercept_calibration(
  batches = 20,
  simulations_per_batch = 500,
  permutations = 99,
  seed = 20260719
)

dispersion_validation <- simulate_dispersion_validation(
  simulations = 500,
  permutations = 99,
  type = "median",
  seed = 20260718
)

hybrid_validation <- simulate_hybrid_validation(
  simulations = 500,
  permutations = 99,
  regions = 8,
  sites = 3,
  replicates = 2,
  seed = 20260718
)

unbalanced_validation <- simulate_unbalanced_validation(
  simulations = 500,
  permutations = 99,
  effect = 0.8,
  seed = 20260718
)

strong_pairwise_validation <- simulate_strong_pairwise_validation(
  simulations = 500,
  permutations = 99,
  effect = 1.2,
  seed = 20260726
)

benchmark_permanova(
  group_counts = c(8, 16, 32),
  permutations = c(99, 999)
)
```

Simulation output includes Monte Carlo uncertainty, denominator agreement,
failure counts, dispersion sensitivity, and runtime. Benchmark results are
machine-specific measurements, not universal performance claims. See the
package vignette for supported-design boundaries and interpretation.

## Vignette series

The guide begins with ordinary ecological questions and introduces the
implementation detail only when random effects make it useful.

- **Getting started with distmer** (`vignette("getting-started")`): what a
  distance-based mixed model is, how it relates to PERMANOVA, and a first
  `vegan::dune` analysis;
- **Fixed effects: Type I, II, and III** (`vignette("fixed-effects")`): choosing
  an estimable fixed-effect question with the `vegan::mite` data;
- **Pairwise contrasts and multiple comparisons** (`vignette("contrasts")`):
  planned, pairwise, marginal, and covariate-adjusted comparisons using `dune`;
- **Dispersion and distance diagnostics** (`vignette("dispersion")`): distance
  geometry, design-aware dispersion, and `vegan::betadisper()` agreement;
- **Unbalanced designs** (`vignette("unbalanced-mixed")`): unequal fixed cells,
  missed repeated visits, unit weighting, and restricted permutations;
- **Random effects and ecological sampling designs**
  (`vignette("random-effects")`): repeated sites, nesting, blocks, hybrids,
  random slopes, and the specialised implementation reference;
- **Validation and reproducibility** (`vignette("validation")`): reference
  agreement, deterministic checks, simulations, and a reporting checklist.

## Example

```r
library(distmer)

Y <- matrix(
  c(
    0, 0,
    0, 1,
    1, 0,
    4, 4,
    4, 5,
    5, 4
  ),
  ncol = 2,
  byrow = TRUE
)

group <- factor(rep(c("A", "B"), each = 3))
D <- dist(Y)

permanova_oneway(
  distance = D,
  group = group,
  permutations = 9999,
  seed = 5,
  engine = "pairwise"
)
```

## Multi-factor development engine

```r
dat <- expand.grid(
  treatment = factor(c("control", "impact")),
  block = factor(seq_len(4)),
  replicate = seq_len(3)
)

set.seed(1)
Y <- matrix(rpois(nrow(dat) * 8, lambda = 4), ncol = 8)
D <- dist(Y)

permanova_fixed(
  distance = D,
  formula = ~ treatment * block,
  data = dat,
  permutations = 999,
  seed = 5
)
```

The `permanova_fixed()` function tests all terms over residual and therefore
must not itself be used for random, nested, repeated, or mixed designs. Use the
dedicated mixed engines only within their documented support boundaries.
