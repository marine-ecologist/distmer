# Validation protocol

Run `run_permanova_validation()` for the deterministic gate,
`simulate_mixed_validation()` for whole-plot operating characteristics, and
`simulate_crossed_interaction_validation()` for the critical
fixed-by-random-interaction calibration experiment,
`simulate_nested_validation()` for top-level balanced `region/site` designs,
and `simulate_nested_level_validation()` for region-, site-, and
observation-level hypotheses. Use
`simulate_crossed_permutation_validation()` for the canonical randomized-block
empirical p-value path and its dispersion-sensitivity scenario.
Use `simulate_dispersion_validation()` for the separate design-aware PERMDISP
companion across all four supported balanced structures.
Use `simulate_hybrid_validation()` for the balanced treatment × region/site
crossed-nested profile-permutation engine.
Use `simulate_unbalanced_validation()` for unequal-cell Type II fixed-effect
inference and its direct response-level Freedman--Lane reference.

The deterministic gate covers:

- Euclidean agreement with ordinary one-way ANOVA;
- balanced random-intercept agreement with ANOVA on whole-plot means;
- exhaustive enumeration of a small whole-plot permutation space;
- preservation and reporting of negative eigenvalues;
- explicit rejection of unsupported imbalance;
- numerical idempotence under near-collinearity; and
- preservation of explicit `:` interaction provenance versus `/` nesting;
- exact agreement between a declared treatment-by-block denominator and the
  classical interaction test;
- rejection of omitted fixed-by-random covariance components by default;
- suppression of Satterthwaite p-values for non-Euclidean spectra; and
- exhaustive within-block crossed-cell permutation and direct statistic
  agreement;
- retained non-Euclidean crossed permutation p-values;
- exact balanced nested agreement with ANOVA on top-unit means;
- exact balanced nested agreement at site and observation levels;
- exhaustive permutation globally, within parents, and within deepest units;
- explicit rejection of fixed terms spanning hierarchy strata; and
- agreement with `vegan::adonis2()` for the supported fixed one-way case when
  `vegan` is installed;
- agreement of dispersion distances with `vegan::betadisper()`;
- exact ANOVA agreement for distances to group centers; and
- exhaustive crossed full-model-residual dispersion statistics against a
  direct blocked-ANOVA reference.

Simulation summaries report null rejection, power, sensitivity to dispersion
heterogeneity, Monte Carlo standard error, denominator agreement, failures,
runtime, and non-finite results. A small number of simulations is not evidence
of nominal Type I error or adequate power.

The crossed-interaction report includes binomial 95% confidence intervals and
compares three results: the declared-interaction engine, the classical
interaction-denominator test, and the explicitly opted-in but misspecified
additive model. Only the first two are expected to agree under generated
interaction variance.

The release-record calibration is pinned and reproducible:

```r
simulate_crossed_interaction_validation(
  simulations = 500,
  blocks = 8,
  replicates = 3,
  interaction_sd = 2,
  block_sd = 1,
  residual_sd = 1,
  alpha = 0.05,
  seed = 20260718
)
```

Its rejection rates are a seeded Monte Carlo realization and will vary with
the seed. The non-stochastic regression requirement is numerical agreement of
the declared-interaction and classical F statistics and p-values.

The nested calibration is likewise pinned and reproducible:

```r
simulate_nested_validation(
  simulations = 500,
  regions_per_treatment = 4,
  sites_per_region = 2,
  replicates = 3,
  region_sd = 1.5,
  site_sd = 1,
  residual_sd = 0.5,
  alpha = 0.05,
  seed = 20260718
)
```

This run gives identical 4.4% rejection for the nested engine and classical
top-unit ANOVA (95% binomial CI 2.8--6.6%). Their F statistics differ by at
most `7.2e-14` and reference p-values by at most `2.4e-13`. The numerical
agreement, rather than one seed's rejection rate, is the regression invariant.

The level-specific calibration is pinned by:

```r
simulate_nested_level_validation(
  simulations = 500,
  regions_per_treatment = 4,
  replicates = 3,
  region_sd = 1.5,
  site_sd = 1,
  residual_sd = 0.5,
  alpha = 0.05,
  seed = 20260718
)
```

Nested and classical tests reject identically at every level: 4.4% for region
(95% CI 2.8--6.6%), 5.0% for site (3.3--7.3%), and 4.2% for observation
(2.6--6.3%). Maximum F differences are `8.9e-14`, `5.9e-14`, and `1.1e-13`;
maximum reference-p differences are `2.9e-13`, `1.9e-13`, and `1.7e-14`.

The crossed randomized-block calibration is pinned by:

```r
simulate_crossed_permutation_validation(
  simulations = 500,
  permutations = 99,
  blocks = 8,
  replicates = 3,
  interaction_sd = 2,
  block_sd = 1,
  residual_sd = 1,
  dispersion_ratio = 2,
  alpha = 0.05,
  seed = 20260718
)
```

Under the exchangeable null, the engine and direct response-level permutation
reference both reject 4.4% (95% CI 2.8--6.6%); the classical interaction-F
reference rejects 5.6% (3.8--8.0%). Engine/direct p-values are identical, and
all permutation statistics agree within `6.9e-13`. With treatment-specific
residual SD ratio 2, the empirical paths reject 5.6% (3.8--8.0%) and the
classical reference 6.6% (4.6--9.1%). This second scenario is a sensitivity
result under violated exchangeability, not evidence of general robustness.

The design-aware dispersion calibration is pinned by:

```r
simulate_dispersion_validation(
  simulations = 500,
  permutations = 99,
  alpha = 0.05,
  seed = 20260718,
  type = "median"
)
```

For the exchangeable null, rejection is 4.6% for one-way designs (95% CI
2.9--6.8%), 6.6% for random-intercept designs (4.6--9.1%), 2.6% for nested
designs (1.4--4.4%), and 5.2% for crossed randomized blocks (3.4--7.5%). A
pure location shift with equal dispersion rejects 4.4%, 4.0%, 3.2%, and 6.0%,
respectively. With residual SD ratio 2, power is 97.4%, 76.6%, 80.6%, and
94.6%. Every observed dispersion F agrees with its independent scalar ANOVA
reference within `8.5e-12`. Centroid simulations are retained as a sensitivity:
they are more liberal in the small one-way and whole-plot designs and are not
the primary calibration claim.

The hybrid calibration is pinned by:

```r
simulate_hybrid_validation(
  simulations = 500,
  permutations = 99,
  regions = 8,
  sites = 3,
  replicates = 2,
  region_sd = 1,
  site_sd = 0.7,
  interaction_sd = 1.2,
  nested_interaction_sd = 0.6,
  residual_sd = 0.4,
  effect = 1,
  alpha = 0.05,
  seed = 20260718
)
```

Under the exchangeable null, engine and direct profile permutation both reject
4.4% (95% CI 2.8--6.6%); the classical interaction-F reference rejects 5.4%
(3.6--7.8%). With fixed location shift 1, permutation power is 24.8%
(21.1--28.8%) and classical power is 27.4% (23.5--31.5%). Engine/direct
p-values are identical. Observed F and every retained permutation statistic
agree within `8.0e-13` and `1.2e-12`, respectively.

The unbalanced fixed-effect calibration is pinned by:

```r
simulate_unbalanced_validation(
  simulations = 500,
  permutations = 99,
  effect = 0.8,
  residual_sd = 1,
  alpha = 0.05,
  seed = 20260718
)
```

For unequal cell counts `3, 5, 4, 2, 6, 3`, the Type II engine and independent
Freedman--Lane reference both reject 5.0% under the null (95% CI 3.3--7.3%);
the classical F reference rejects 4.0% (2.5--6.1%). With fixed A effect 0.8,
all three paths reject 41.0% (36.7--45.5%). Engine/direct p-values are
identical. Maximum observed-F and permutation-statistic differences are
`6.4e-14` and `9.6e-14`.

The unequal-replication random-intercept calibration is pinned by:

```r
simulate_unequal_random_intercept_validation(
  simulations = 500,
  permutations = 99,
  effect = 1,
  cluster_sd = 1,
  residual_sd = 1,
  alpha = 0.05,
  seed = 20260719
)
```

For unit sizes 2, 3, and 5, each represented by two units per treatment, the
cluster-centroid engine and independent scalar reconstruction reject 7.6%
under the pinned null (95% CI 5.4--10.3%) and have 22.6% power for shift 1
(19.0--26.5%). The mildly high seeded null realization is reported as observed,
not as a universal calibration constant. The regression invariant is exact
engine/direct agreement: p-values are identical, observed F differs by at most
`2.8e-13`, and permutation statistics by at most `2.3e-13`.

The unbalanced nested/randomized-block calibration is pinned by:

```r
simulate_unbalanced_mixed_validation(
  simulations = 500,
  permutations = 99,
  effect = 1,
  parent_sd = 1,
  child_sd = 0.7,
  block_sd = 1,
  interaction_sd = 1,
  residual_sd = 0.5,
  alpha = 0.05,
  seed = 20260719
)
```

Nested engine/direct null rejection is 3.4% (95% CI 2.0--5.4%) and location
power is 24.4% (20.7--28.4%). Incomplete-block engine/direct null rejection is
6.4% (4.4--8.9%) and power is 28.4% (24.5--32.6%). Engine/direct p-values are
identical across both structures; maximum observed-F and permutation-statistic
differences are `1.2e-13` and `2.2e-13`.

The planned-contrast calibration is pinned by:

```r
simulate_contrast_validation(
  simulations = 500,
  permutations = 99,
  effect = 0.8,
  residual_sd = 1,
  alpha = 0.05,
  seed = 20260719
)
```

For unequal group counts 4, 6, and 5, the planned A-minus-B engine/direct null
rejection is 5.6% (95% CI 3.8--8.0%) and the classical one-df reference rejects
5.2% (3.4--7.5%). Permutation power for shift 0.8 is 18.6% (15.3--22.3%); the
classical reference has 20.8% power (17.3--24.6%). Engine/direct p-values are
identical. Maximum engine/direct observed-F and permutation-statistic
differences are `5.7e-14` and `2.9e-13`; observed F agrees with the classical
contrast formula within `5.7e-14`.

The Phase 14 unequal-random-intercept calibration audit is pinned by:

```r
audit_unequal_random_intercept_calibration(
  batches = 20,
  simulations_per_batch = 500,
  permutations = 99,
  engine_checks = 10,
  cluster_sd = 1,
  residual_sd = 1,
  alpha = 0.05,
  seed = 20260719
)
```

Each scenario contains 10,000 datasets. Gaussian balanced allocation rejects
4.76% exactly (95% CI 4.35--5.20%) and 4.83% with 99 Monte Carlo draws
(4.42--5.27%). Gaussian unequal allocation rejects 4.58%/4.95%; heavy-tailed
random effects and residuals reject 4.53%/4.42%. The treatment-specific
dispersion scenario has `exchangeability = FALSE` and does not contribute to
the nominal calibration gate.

The first Gaussian batch remains high at 6.8% exact and 6.6% under the audit's
Monte Carlo draws, while the median exact rejection across 20 independent
batches is 4.8%. This demonstrates why the earlier full-engine 7.6% batch must
not be promoted to a general operating characteristic. Batch seeds are spaced
one million apart; adjacent base seeds would reuse nearly all deterministic
per-iteration streams and are not independent evidence.

The full distance engine matches the unique-label assignment reference within
`8.9e-16` for F and exactly for p-values. In the balanced-allocation design,
`(4!)^3 = 13,824` distinct profile rearrangements collapse to
`choose(4, 2)^3 = 216` unique label assignments. The intercept-only
Freedman--Lane transformation is algebraically identical to raw centroid
permutation, and this identity is part of the deterministic gate.

The Phase 15 marginal-contrast validation is pinned by:

```r
simulate_marginal_contrast_validation(
  simulations = 500,
  permutations = 99,
  effect = 0.8,
  residual_sd = 1,
  alpha = 0.05,
  seed = 20260719
)
```

In the unequal-cell `A * B` design, the equal-grid marginal A contrast rejects
4.8% under the null (95% CI 3.1--7.1%) for both the engine and independent
response-level reconstruction; the classical GLH rejects 4.0% (2.5--6.1%).
For shift 0.8, permutation power is 38.8% (34.5--43.2%) and classical power is
40.0% (35.7--44.4%). Engine/direct p-values are identical. Maximum observed-F
and permutation-statistic differences are `4.3e-14` and `2.2e-13`.

The Phase 16 mixed-contrast validation is pinned by:

```r
simulate_mixed_contrast_validation(
  simulations = 500,
  permutations = 99,
  blocks = 8,
  replicates = 2,
  effect = 0.8,
  alpha = 0.05,
  seed = 20260719
)
```

For the declared three-treatment randomized-block model, the planned A-minus-B
contrast rejects 5.0% under the null by empirical inference (95% CI
3.3--7.3%); the Gaussian scalar GLH rejects 5.4% (3.6--7.8%). For shift 0.8,
empirical power is 26.6% (22.8--30.7%) and scalar-reference power is 27.8%
(23.9--32.0%). The observed mixed-contrast F equals the independent scalar GLH
within `1.2e-13`. Deterministic tests additionally cover all four structures
and exhaustive response-level permutation-vector equality.

The Phase 17 global-null max-F validation is pinned by:

```r
simulate_maxT_validation(
  simulations = 5000,
  permutations = 99,
  effect = 0.8,
  alpha = 0.05,
  seed = 20260719
)
```

Across 5,000 unequal-size three-group null datasets, unadjusted pairwise family
rejection is 11.36% (95% CI 10.49--12.27%), Holm is 2.68%
(2.25--3.17%), and max-F is 4.92% (4.34--5.56%). Under a group-C shift of 0.8,
family rejection is 34.00%, 11.72%, and 19.48%, respectively. The invariant
validation is exact brute-force equality of the complete joint statistic and
maximum distributions; the operating-characteristic rate checks only weak
FWER under the complete null.

## Stratified-span safeguard amendment

The Phase 11/12 safety amendment pins a numerical boundary that relative QR
rank cannot detect reliably. If treatment is constant within every cluster-
size, hierarchy-topology, or block-by-cell-size stratum, stratum centering
should remove its entire estimable direction. Floating arithmetic may leave a
projected singular value around `1e-16`; a relative QR tolerance can rescale
that dust and report rank one.

The amended helper compares projected singular values with an absolute floor
defined relative to the unprojected treatment-design scale. Adversarial tests
now require explicit rejection for fully confounded examples in all three
engines. Valid Phase 11 and Phase 12 F statistics and exhaustive permutation
vectors remain equal to their independent unit-centroid references.

## Phase 19 unbalanced mixed contrasts

The pinned validation command is:

```r
simulate_unbalanced_mixed_contrast_validation(
  simulations = 500,
  permutations = 99,
  effect = 0.8,
  alpha = 0.05,
  seed = 20260720
)
```

The unequal-size three-treatment random-intercept design has every treatment
represented within exact cluster sizes 2, 3, and 4. Engine and direct
restricted-permutation null rejection are both 5.0% (95% CI 3.3--7.3%); the
Gaussian scalar reference rejects 2.8% (1.5--4.7%). Under shift 0.8,
permutation power is 12.0% (9.3--15.2%) and scalar-reference power is 11.6%
(8.9--14.7%). Engine/direct p-values are identical. Maximum observed-F and
permutation-vector differences are `2.3e-13` and `2.8e-12`.

Deterministic references separately cover an unbalanced parent/child hierarchy
and an incomplete randomized block. The latter averages fitted treatment means
over one common equal-weight block grid, preventing observed block incidence
from leaking into the treatment contrast.

## Phase 20 covariate-adjusted contrasts

The pinned validation command is:

```r
simulate_adjusted_contrast_validation(
  simulations = 500,
  permutations = 99,
  effect = 0.8,
  alpha = 0.05,
  at = 0.25,
  seed = 20260721
)
```

The unbalanced three-treatment fixed design contains one categorical nuisance
factor, one continuous covariate, and a treatment-by-covariate interaction.
Engine and direct Freedman--Lane null rejection are both 4.0% (95% CI
2.5--6.1%); the Gaussian scalar GLH rejects 4.2% (2.6--6.3%). Under shift 0.8,
permutation power is 11.2% (8.6--14.3%) and scalar-reference power is 10.8%
(8.2--13.9%). Engine/direct p-values are identical. Maximum observed-F and
permutation-vector differences are `7.4e-13` and `2.0e-13`.

The invariant checks are coefficient-row equality with an explicit reference
grid, observed-F equality with a scalar GLH, and direct reconstruction of the
sampled residual-permutation statistics. The operating-characteristic rates
are a reproducible audit, not fixed package constants.

## Phase 21 mixed-design global-null max-F

The pinned weak-FWER audit is:

```r
simulate_mixed_maxT_validation(
  simulations = 5000,
  permutations = 99,
  effect = 0.8,
  alpha = 0.05,
  seed = 20260722
)
```

The unequal-size three-treatment random-intercept design represents every
treatment within exact cluster sizes 2, 3, and 4. Under the complete null,
unadjusted pairwise family rejection is 13.36% (95% CI 12.43--14.33%), Holm is
1.84% (1.49--2.25%), and max-F is 3.82% (3.31--4.39%). Under a group-C shift of
0.8, family rejection is 21.00%, 3.36%, and 6.38%, respectively.

Observed statistics agree with independently projected scalar unit-centroid
GLHs within `5.2e-13`; all sampled joint statistic matrices agree within
`3.5e-12`; raw and maximum p-values are identical. The max-F rate audits weak
FWER only under the complete treatment null and does not establish strong
subset-null control.

## Phase 22 v2 consolidation gate

`permanova_capabilities()` is checked for 17 unique primary inference engines,
resolvable pairwise companions, complete design/estimand/denominator/
permutation/error-control/exclusion fields, and explicit weak-FWER labels on
both maximum-statistic engines. `inst/SUPPORT-MATRIX.md` must name every
registered engine and preserve the statement that runtime validators are
authoritative.

The final v2 release gate runs, in order:

```sh
Rscript tools/release-audit.R
Rscript -e 'testthat::test_local(reporter = "progress")'
R CMD build .
R CMD check --no-manual distmer_2.0.0.tar.gz
R CMD check --as-cran --no-manual distmer_2.0.0.tar.gz
```

The offline `--as-cran` triage may retain environment-only CRAN connectivity or
clock notes. Any package warning or error remains a failed gate. Owner identity,
copyright, URL, and bug-tracker values are publication blockers and are not
invented by the implementation.

The completed v2.0.0 gate passes 57/57 deterministic checks and 682 package
expectations. All 20 vignettes build. The normal check reports `Status: OK`;
the offline `--as-cran` check reports one system-clock NOTE (`unable to verify
current time`) and no package warning or error. Publication metadata remains
the only intentional release blocker.

## Phase 23 external adversarial validation

Claude independently verified the pinned SHA-256 values for the sealed v2.0.0
and Phase 25 development artifacts, then attacked Phases 19--25. The adjusted
two-level contrast pseudo-F agreed with scalar ANCOVA within `4.8e-14`.
Mixed max-F recomputation matched the declared permutation algebra within
`1e-12`; in the independent complete-null audit its weak family-error rate was
5.0% (95% CI 3.8--6.4%) versus 11.8% for unadjusted minimum p-values. The lower
unbalanced-nested child-level statistic agreed with its equal-weight centroid
reference within `1.1e-14`, and the parent-by-size confounding attack rejected
through the shared absolute-scale span guard.

One diagnostic gap was reproduced: a numerically full-rank but effectively
aliased covariate could reach `qr.solve(crossprod(X), ...)` and expose a LAPACK
singularity message. The development line now checks the smallest-to-largest
singular-value ratio of `X` against `sqrt(tol)` before projector construction.
The exact adversarial design rejects as ill-conditioned, while a correlated
but stable control remains supported. This check is part of the deterministic
gate and the ratio is returned by `derive_adjusted_means()`.

The amended development gate passes 62/62 deterministic checks and 722/722
test expectations; all 21 vignettes build and the normal package check reports
`Status: OK`. Offline `--as-cran` reports only the unavailable incoming check,
the intentional development-version NOTE, and the environment clock NOTE. The
amended artifact is
`/tmp/distmer_2.1.0.9000-phase23-amended.tar.gz`, SHA-256
`e5a151e9b7533af30189e6c2838f9a7c8f6456243395e74c4915b8d77b0abd54`.

## Phase 25 lower-level unbalanced nested inference

The pinned validation command is:

```r
simulate_unbalanced_nested_lower_validation(
  simulations = 500,
  permutations = 99,
  effect = 0.8,
  alpha = 0.05,
  seed = 20260723
)
```

The child-level design uses equally weighted child centroids, parent
adjustment, and parent-by-exact-size residual permutations. Engine and direct
null rejection are both 6.0% (95% CI 4.1--8.5%); permutation power under shift
0.8 is 39.2% (34.9--43.6%). The observation-level design adjusts for child and
permutes residual observations within child. Engine/direct null rejection is
3.8% (2.3--5.9%) and power is 95.2% (92.9--96.9%).

Engine/direct p-values are identical at both levels. Maximum observed-F and
permutation-vector differences from the scalar/direct references are
`2.5e-13` and `1.4e-13`. Deterministic tests additionally cover exhaustive
enumeration, lower-level planned/pairwise contrasts, mixed max-F routing,
wrong requested levels, disconnected incidence, parent-by-size confounding,
zero residual df, zero denominators, and retained non-Euclidean inference.

The completed development gate passes 61/61 deterministic checks and 719
package expectations. All 21 vignettes build, the normal check reports
`Status: OK`, and offline `--as-cran` reports only the intentional
2.1.0.9000-version and environment clock NOTEs, with no package warning or
error.

## Phase 26 unbalanced crossed-nested hybrid inference

The pinned validation command is:

```r
simulate_unbalanced_hybrid_validation(
  simulations = 500,
  permutations = 99,
  effect = 0.8,
  alpha = 0.05,
  seed = 20260724
)
```

The audit uses unequal sites per region, unequal replicated treatment-site
cells, and one missing whole treatment-region profile. Every observed profile
contains its region's complete site set and averages site-cell centroids with
equal weight. Engine and direct null rejection are both 3.8% (95% CI
2.3--5.9%); scalar parametric rejection is 3.6% (2.1--5.6%). Under shift 0.8,
engine/direct power is 14.4% (11.4--17.8%) and scalar power is 16.6%
(13.4--20.2%).

Engine and direct p-values are identical. Across 1,000 generated datasets,
maximum pseudo-F disagreement from the scalar profile GLH is `4.91e-13`, and
the maximum sampled permutation-vector difference is `5.63e-12`. Deterministic
tests additionally cover exhaustive 432-arrangement reconstruction, common-
region mixed contrasts, shared max-F routing, balanced-EMS routing, partial
profile rejection, disconnected incidence, topology confounding, a zero
profile residual, and retained non-Euclidean inference.

The completed `2.2.0.9000` gate passes 64/64 deterministic checks and 749/749
test expectations. All 22 vignettes build and the normal check reports
`Status: OK`; offline `--as-cran` has only the unavailable-incoming/
development-version and environment-clock NOTEs. The pinned artifact is
`/tmp/distmer_2.2.0.9000-phase26.tar.gz`, SHA-256
`ee5d2efc49a163f15d49021a3247961538e27b9937bcebfbaeb48930662d2d19`.

## Phase 27 correlated random intercepts and slopes

The pinned validation command is:

```r
simulate_random_slope_validation(
  simulations = 500,
  permutations = 99,
  effect = 0.8,
  alpha = 0.05,
  seed = 20260725
)
```

Under the centrally symmetric null, engine, direct sign-flip, and Gaussian
paired-difference references each reject 4.4% (95% CI 2.8--6.6%). Under shift
0.8, engine/direct power is 56.4% (51.9--60.8%) and the Gaussian reference is
63.6% (59.2--67.8%). Engine/direct p-values are identical. Maximum scalar-F,
sign-flip-vector, and covariance-matrix differences are respectively
`1.71e-13`, `1.03e-12`, and `3.11e-15`.

A centred asymmetric random-slope distribution is reported only as a failed-
assumption sensitivity: sign-flip rejection is 7.4% (5.3--10.1%), whose lower
confidence bound exceeds 5%. This is positive evidence for the boundary, not
an extension of the guarantee. Deterministic tests also cover exact `2^n`
enumeration, paired-t equality, residual-corrected covariance, reproducible
Monte Carlo signs, multi-level/missing/unequal/unreplicated rejection, zero
denominators, and retained semimetric inference.

The completed `2.3.0.9000` gate passes 67/67 deterministic checks and 786/786
test expectations. All 23 vignettes build; the normal installed-package check
reports `Status: OK`, and the offline `--as-cran` check has no package warning
or error. The pinned artifact is
`/tmp/distmer_2.3.0.9000-phase27.tar.gz`, SHA-256
`a03c6190423277baa0728d25dc89645ea3e5ea5f9b30fa4cd80b15bd9ac7a33a`.

## Phase 28 strong closed pairwise inference

The existing shared max-F engines remain weak complete-null procedures.
Phase 28 proceeds only for a family with a direct subset-null invariance proof:
all pairwise complete-distribution equalities among levels of one fixed factor.
For every intersection, the asserted equality edges define connected
components; raw observations are reassigned only inside those components with
group sizes preserved, and all outside groups remain fixed. Equality of the
complete distributions inside each component makes this action invariant.
Local maximum-F tests combined by closure therefore provide strong FWER
without a subset-pivotality claim.

The pinned validation command is:

```r
simulate_strong_pairwise_validation(
  simulations = 500,
  permutations = 99,
  effect = 1.2,
  alpha = 0.05,
  seed = 20260726
)
```

Complete-null strong FWER is 3.4% (95% CI 2.0--5.4%). With `A` and `B`
identically distributed and `C` shifted by 1.2, rejection of the true `A = B`
pair is 2.0% (1.0--3.6%) and power against at least one false pair is 42.0%.
Across all 1,000 datasets, every sampled intersection statistic and every
local/closed p-value agrees exactly with direct reconstruction, and the
partial-null `A = B` action moves zero `C` observations.

The guarantee is not centroid-only validity under heterogeneous distributions.
Planned contrasts, dependent observations, repeated measures, and mixed
designs remain outside the Phase 28 proof and reject rather than inheriting a
strong-control label.

The completed `2.4.0.9000` gate passes 70/70 deterministic checks and 835/835
test expectations. All 24 vignettes build; the normal installed-package check
reports `Status: OK`, and offline `--as-cran` has no package warning or error.
The pinned artifact is `/tmp/distmer_2.4.0.9000-phase28.tar.gz`,
SHA-256
`ad8aed00d13f43859ea4130ee1c67137dcd9c51b2f3a94b9cd6d35e2ad04c5ce`.

## Pre-1.0 mixed-inference adversarial hardening

The post-Phase-28 development review added regression anchors for every defect
reproduced in the unified, covariate, Monte-Carlo, and continuous-slope paths:

- a between-unit numeric predictor is assigned to the between-group stratum in
  a single-random-intercept model, while the nested covariate engine rejects it
  because no supported deepest-level numerator exists;
- Type I and Type III balanced mixed tests reproduce their distinct classical
  whole-plot references in a non-orthogonal factorial design;
- nested covariate Type III uses the interaction-adjusted partial projector,
  and the Freedman--Lane reduced projector contains the same nuisance terms;
- exact mixed-covariate actions enumerate 24 whole-group arrangements and 16
  within-group arrangements in the canonical four-group/two-observation case;
- scaling only zero-mean within-site deviations leaves both an upper-stratum
  pseudo-F and its `P(MC)` unchanged; and
- continuous random-slope parametric inference is suppressed outside rank-one,
  equal-`SSx` geometry, while its sampling correction uses
  `MS[Residual] * mean(1 / SSx)`.

The completed `0.9.0.9000` gate passes 70/70 deterministic checks and 932/932
source test expectations with no failures, errors, warnings, or skips. All seven
ecologist-oriented vignettes build, and `R CMD check --no-manual` reports
`Status: OK`.
The reproducible source artifact is built as
`/private/tmp/distmer_0.9.0.9000.tar.gz`; its checksum is recorded externally
after the final build so the validation document does not contain a
self-invalidating hash of its own artifact.
At that 0.9.0.9000 checkpoint, publication metadata remained unresolved. Those
fields are present and validated in the current gate.

## Version 0.9.2 consolidated Claude development and adversarial gate

Version 0.9.2 incorporates the complete post-0.9.1 Claude development tree:
the shared low-rank Freedman--Lane kernel, asymmetric three-stratum mixed
inference, balanced and unbalanced multi-factor split-plots, fit-based
contrasts, progress reporting, and their generated documentation.

The independent adversarial gate adds direct anchors for 24 exact whole-unit
and 576 exact within-unit asymmetric assignments; incomplete units restricted
by subplot-count topology; 24 whole-unit and 16 within-unit split-plot
assignments; exact-size conditioning for unequal unit centroids; and mixed
isolated contrasts reproducing the complete subset model.

The checked source passes 69/69 deterministic checks and 1079/1079 expectations
across 239 tests, with no failures or skips. Both release audits pass, all seven
vignettes build, and `R CMD check --no-manual` reports `Status: OK`.
