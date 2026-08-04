# distmer: Transparent Distance-Based Models for Ecological Communities

Transparent distance algebra, canonical design compilation, and
permutation inference for ecological community models. Initial scope is
one-way and sequential fixed-effect PERMANOVA, canonical mixed-design
compilation, and exact expected-mean-square inference with restricted
permutations for balanced single-random-intercept models. Balanced
crossed models add trace-derived EMS, declared fixed-by-random
interaction denominators, quasi-denominators, and unrestricted
distance-based variance-component estimates. Omitted interaction
components reject unless additivity is explicitly requested. The
canonical balanced fixed-treatment by random-block design also supports
within-block permutation of intact replicated treatment cells, including
for non-Euclidean dissimilarities. Balanced nested random-intercept
hierarchies add exact top-unit denominators, intact-hierarchy residual
permutations, and trace-derived nested variance components. Fixed terms
may be assigned to top, child, or observation strata when each estimable
hypothesis lies wholly in one balanced hierarchy level. Design-aware
multivariate dispersion tests use spatial medians or centroids, retain
the full pseudo-Euclidean spectrum, and reuse the validated balanced EMS
and restricted residual-permutation structures. Balanced crossed-nested
hybrids support one treatment crossed with regions and sites nested
within regions, with explicit treatment-by-hierarchy covariance
components and intact treatment-profile permutations. Unbalanced
fixed-effect models add explicit estimable Type I, II, or III
hypotheses, alias diagnostics, full-model residual denominators, and
term-specific residual permutations. Unequally replicated
single-random-intercept designs use equally weighted cluster centroids
and whole-cluster residual permutations restricted by exact cluster
size. Unbalanced parent/child hierarchies add topology-stratified
top-centroid inference, while incomplete randomized blocks add intact
cell-centroid permutations restricted jointly by block and cell size.
Unbalanced crossed-nested region/site hybrids use equal-site-weight
treatment-region profile centroids, additive region adjustment, and
empirical permutations restricted by exact site-specific replication
topology; no unbalanced EMS approximation is made. Balanced complete
two-condition repeated-measures designs add correlated random
intercept–slope trace covariance estimates and whole-subject sign-flip
inference under an explicit central-symmetry assumption. Child-level
treatments in unbalanced parent/child hierarchies use equally weighted
child centroids, parent adjustment, and parent-by-size restricted
residual permutations. Observation-level treatments adjust for child and
permute residual observations within child. Planned and pairwise
contrasts among one-factor centroids use contrast-specific residual
permutations and explicit familywise or false-discovery-rate adjustment.
Exact unique-assignment and non-overlapping multi-batch audits harden
unequal-random-intercept calibration under Gaussian, heavy-tailed, and
unequal-allocation nulls. Full-rank categorical fixed models add equal-
or frequency-weighted marginal means, planned contrasts, and adjusted
pairwise comparisons. Full-rank fixed models with ordinary numeric
covariates add contrasts at explicit covariate reference values and
declared categorical nuisance weights. Balanced mixed designs add
planned and pairwise contrasts on the correct experimental-unit centroid
scale, with whole-unit or within-parent restricted residual permutations
and exact interaction denominators. Unequally replicated
random-intercept, unbalanced top-level nested, and incomplete
randomized-block designs extend those contrasts to explicitly stratified
empirical unit-centroid inference. Fixed one-factor contrast families
add shared global-null max-F permutation adjustment with explicitly weak
familywise error control. Supported mixed designs extend that weak
complete-null max-F adjustment to shared raw experimental-unit
permutations within their declared size, topology, block, or region
strata. Fixed one-factor all-pair families add strong familywise control
for complete-distribution equality nulls by closed testing with a
separately valid connected-component label-permutation action at every
intersection. This guarantee does not extend to centroid-only nulls
under heterogeneous distributions. A machine-readable capability
registry consolidates every validated engine's estimand, denominator,
permutation contract, error-control claim, and explicit exclusions.

## See also

Useful links:

- <https://marine-ecologist.github.io/distmer/>

- <https://github.com/marine-ecologist/distmer>

- Report bugs at <https://github.com/marine-ecologist/distmer/issues>

## Author

**Maintainer**: Roff George <groff@reefresearch.org>

Authors:

- Roff George <groff@reefresearch.org>
