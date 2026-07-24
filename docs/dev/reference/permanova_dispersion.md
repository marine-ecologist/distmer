# Design-aware multivariate dispersion inference

Converts the original dissimilarity to distances from each observation
to its target-group spatial median or centroid, then tests that scalar
response with the EMS denominator and restricted residual permutation
appropriate to a validated balanced design. The initial envelope
contains exactly one categorical fixed main effect.

## Usage

``` r
permanova_dispersion(
  formula = NULL,
  data = NULL,
  distance = NULL,
  term = NULL,
  type = c("median", "centroid"),
  bias_adjust = FALSE,
  permutations = 999L,
  seed = NULL,
  keep_permutations = FALSE,
  max_exact = 100000L,
  tol = 1e-10,
  design = NULL
)
```

## Arguments

- formula:

  A model formula, or a compiled \`permanova_design\` object. The
  dissimilarity may be supplied on the left-hand side, for example
  \`bray ~ treatment + (1 \| site)\`, or through \`distance\` with a
  one-sided formula. A compiled design remains accepted for low-level
  workflows.

- data:

  A data frame with one row per observation. Do not supply it when
  \`formula\` is already a compiled design. For backward compatibility,
  the second positional argument after a compiled design may be the
  target term name.

- distance:

  A \`dist\` object or full distance matrix, used with a one-sided
  formula.

- term:

  Optional name of the sole categorical fixed main-effect term. It is
  inferred when the model contains exactly one fixed term.

- type:

  Group spatial \`"median"\` or \`"centroid"\`.

- bias_adjust:

  Apply the small-sample group-size adjustment.

- permutations:

  A non-negative integer or \`"all"\`.

- seed:

  Optional permutation seed.

- keep_permutations:

  Retain permutation statistics.

- max_exact:

  Maximum exact permutation space.

- tol:

  Numerical tolerance.

- design:

  Backward-compatible named alias for a compiled Design Object; use
  \`formula = design\` or pass the object as the first argument.

## Value

An object of class \`permanova_dispersion\`.

## Examples

``` r
data(dune, dune.env, package = "vegan")
selected <- dune.env$Management %in% c("NM", "SF")
env <- droplevels(dune.env[selected, , drop = FALSE])
bray <- vegan::vegdist(dune[selected, ], method = "bray")
permanova_dispersion(
  bray ~ Management, env, permutations = 19, seed = 1
)
#> Design-aware multivariate dispersion test
#> Term: Management  Center: median 
#> Engine: one-way 
#> Permutation: full-model residual permutation of observations 
#>        term df         ss         ms  pseudo_f p_perm
#>  Management  1 0.02578244 0.02578244 0.9172607   0.25
#>    Residual 10 0.28108080 0.02810808        NA     NA
```
