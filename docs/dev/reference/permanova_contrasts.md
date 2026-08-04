# Contrasts for a fitted PERMANOVA

Planned or pairwise contrasts on a fixed factor of a fitted
[`permanova()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova.md)
model. Contrast families mirror emmeans. Two testing methods are offered
(see `method`), differing in denominator and permutation scope.

## Usage

``` r
permanova_contrasts(
  fit,
  type = "pairwise",
  method = c("isolated", "pooled"),
  factor = NULL,
  ref = NULL,
  scores = NULL,
  custom = NULL,
  adjust = "holm",
  permutations = NULL,
  seed = NULL
)
```

## Arguments

- fit:

  A `permanova` object.

- type:

  Contrast family: "pairwise", "trt.vs.ctrl", "poly", "contr.sum",
  "helmert", or a contrast matrix / named list.

- method:

  Testing method. `"isolated"` (default) runs each contrast as a
  two-group
  [`permanova()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova.md)
  on the involved levels only - its own residual and restricted
  permutation; this matches PRIMER's pairwise routine and reuses the
  validated engine. `"pooled"` tests against the full-model residual
  with numerator SS = trace(H_c G); it borrows error df (more power) but
  assumes homogeneous dispersion, and its permutation null requires the
  engine's own generator (see Details).

- factor:

  Name of the fixed factor. Defaults to the tested factor.

- ref:

  Reference level for "trt.vs.ctrl".

- scores:

  Numeric spacing for "poly".

- custom:

  Matrix / named list for type = "custom".

- adjust:

  Multiplicity method for
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) (default
  "holm").

- permutations:

  Number of permutations; defaults to `fit`'s value.

- seed:

  RNG seed; defaults to `fit`'s value.

## Value

An object of class `permanova_contrasts`.

## Details

`"isolated"` is exact and fully validated: it is literally
[`permanova()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova.md)
on the subset, so its statistic, denominator and permutation match a
two-group analysis of those levels. `"pooled"` uses a numerator SS that
is validated deterministically (it equals the between-group SS to
machine precision), but its permutation null must be produced by the
engine's Freedman-Lane generator restricted to the contrast's
exchangeable units. Until that generator is wired in, `"pooled"` returns
the statistic, df and denominator but leaves the p-value `NA` rather
than computing it over the wrong (full-data) permutation space.

## See also

[`permanova()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova.md),
[`contrast_matrix()`](https://marine-ecologist.github.io/distmer/dev/reference/contrast_matrix.md)
