# Fit a distance-based PERMANOVA for any supported design

A single front door that compiles a formula into a canonical Design
Object (or accepts one directly), routes it to the engine matching its
random-effect provenance and balance, and returns a result that prints
as a PERMANOVA table. No denominator or permutation scheme is ever
chosen silently: a design outside every supported envelope stops with a
message naming the reason.

## Usage

``` r
permanova(
  formula,
  data = NULL,
  distance = NULL,
  type = NULL,
  permutations = 999L,
  mc_permutations = 9999L,
  seed = NULL,
  keep_permutations = FALSE,
  truncate = FALSE,
  max_exact = 100000L,
  tol = 1e-10
)
```

## Arguments

- formula:

  A model formula, or a compiled \`permanova_design\` object. The
  dissimilarity may be given on the left-hand side, as in \`bray ~
  treatment + (1 \| site)\` (the idiom used by \`vegan::adonis2()\`), or
  the formula may be one-sided (\`~ treatment + (1 \| site)\`) with the
  dissimilarity passed through \`distance\`.

- data:

  A data frame with one row per observation. Ignored when \`formula\` is
  already a Design Object.

- distance:

  A \`dist\` object or distance matrix. Needed only when the
  dissimilarity is not given on the formula's left-hand side.

- type:

  Sum-of-squares partition for the fixed-effect terms: \`"I"\`
  (sequential), \`"II"\` (conditional), or \`"III"\` (partial). Applies
  to fixed-effect designs; for designs whose denominator is fixed by the
  random structure it is reported as not applicable.

- permutations:

  A non-negative integer or \`"all"\`.

- mc_permutations:

  Number of Monte-Carlo draws for the asymptotic \`P(MC)\` p-value,
  reported alongside \`P(perm)\` for designs tested in the full
  observation space. \`P(MC)\` stays smooth when the true permutation
  space is small (the usual case for the upper strata of a mixed
  design). Set to \`0\` to disable. Reproducible from \`seed\`.

- seed:

  Optional random seed. Also seeds the \`P(MC)\` draws.

- keep_permutations:

  Retain permutation statistics on the fitted object.

- truncate:

  Also report zero-truncated variance components where the engine
  estimates them.

- max_exact:

  Maximum permutation space allowed for exact enumeration.

- tol:

  Numerical tolerance.

## Value

An object of class \`permanova\` wrapping the engine result. Use
\`summary()\` for the full table with variance components.

## Examples

``` r
data(dune, dune.env, package = "vegan")
bray <- vegan::vegdist(dune, method = "bray")
permanova(bray ~ Management, dune.env, permutations = 99, seed = 1)
#> Distance-based PERMANOVA
#> Engine: fixed effects  (Type I SS)
#> 
#> Factors
#>  Name       Type  Levels
#>  Management Fixed 4     
#> 
#> PERMANOVA table of results
#>  Source     df SS    MS    Pseudo-F P(perm) P(MC)  Unique perms
#>  Management 3  1.469 0.490 2.7672   0.0100  0.0051 99          
#>  Residual   16 2.830 0.177                                     
#>  Total      19 4.299                                           
```
