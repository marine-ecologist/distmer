# Compile a canonical distance-model design object

Parses fixed effects, random-intercept terms, and the narrow correlated
random-slope syntax used by Phase 27 without performing any inference.
The resulting object is intended to be the sole model-structure input to
later EMS and permutation engines.

## Usage

``` r
compile_permanova_design(formula, data, distance = NULL, tol = 1e-10)
```

## Arguments

- formula:

  An \`lme4\`-style formula. Random effects support \`(1 \| group)\`,
  nesting with \`/\`, explicit grouping interactions with \`:\`, and one
  correlated \`(1 + slope \| group)\` declaration.

- data:

  A data frame with one row per observation.

- distance:

  Optional distance object used for row-count and spectral diagnostics.

- tol:

  Numerical rank and eigenvalue tolerance.

## Value

An object of class \`permanova_design\`.
