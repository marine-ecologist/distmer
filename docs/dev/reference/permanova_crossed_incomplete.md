# Crossed mixed PERMANOVA for incomplete, unequally replicated designs

Fixed and random factors fully crossed, where cells may be empty and
replication may be unequal. Sums of squares are marginal: each term is
tested for what it adds to the model containing every other term.
Denominators come from the design's expected mean squares, so a fixed
term crossed with a random factor is tested over the stratum carrying
the same random components, and falls back to a Satterthwaite
combination when imbalance leaves no single stratum matching.

## Usage

``` r
permanova_crossed_incomplete(
  design,
  hypothesis = c("II", "I", "III"),
  permutations = 999L,
  seed = NULL,
  keep_permutations = FALSE,
  max_exact = 100000L,
  tol = 1e-10
)
```

## Arguments

- design:

  A `permanova_design` from
  [`compile_permanova_design()`](https://marine-ecologist.github.io/distmer/dev/reference/compile_permanova_design.md)
  with no `/`-generated nesting. Random terms are optional: with none,
  every term is tested over the residual, which is how an incomplete
  *fixed* design is handled.

- hypothesis:

  Sums of squares. `"II"` (default) adjusts each term for every term
  that does not contain it. `"I"` is sequential, in the model matrix's
  term order. `"III"` adjusts each term for every other term in the
  model; on an incomplete design a term containing this one can span it,
  so Type III may report fewer estimable degrees of freedom than Type II
  – and in a saturated random model can drive a main effect to zero.
  Type II is the default because it stays well posed however incomplete
  the design is.

- permutations:

  Number of Freedman-Lane permutations, or `"all"`.

- seed:

  Optional random seed.

- keep_permutations:

  Retain the permutation statistics.

- max_exact:

  Largest permutation space to enumerate exactly.

- tol:

  Numerical tolerance for ranks.

## Value

An object of class `permanova_crossed_incomplete`.

## Details

Degrees of freedom are estimable dimensions, not level-count products. A
term whose cells are partly empty reports the rank it actually spans,
and `nominal_df` in the returned table records what a complete design
would have given, so the loss is visible rather than silent.

Expected mean squares follow the restricted parameterisation – fixed
effects summing to zero over mixed terms – so a random stratum carries
no component from a random-by-fixed interaction sitting inside it.

## See also

[`permanova_crossed()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_crossed.md)
for complete, equally replicated designs.
