# Balanced distance-based random-intercept inference

Builds an explicit test specification for every fixed term, selects the
exact between-whole-plot residual denominator from its EMS, and applies
hypothesis-specific whole-plot residual permutations.

## Usage

``` r
permanova_mixed(
  design,
  permutations = 999L,
  seed = NULL,
  hypothesis = "I",
  keep_permutations = FALSE,
  max_exact = 100000L,
  tol = 1e-10
)
```

## Arguments

- design:

  A canonical \`permanova_design\` compiled with \`distance\`.

- permutations:

  A non-negative integer or \`"all"\`. Enumeration is used automatically
  when the requested sample covers the attainable space.

- seed:

  Optional random seed.

- hypothesis:

  Fixed-effect sum-of-squares partition: \`"I"\` (sequential), \`"II"\`
  (marginal), or \`"III"\` (partial).

- keep_permutations:

  Retain permutation statistics in each test spec.

- max_exact:

  Maximum permutation space allowed for enumeration.

- tol:

  Numerical tolerance.

## Value

An object of class \`permanova_mixed\`.
