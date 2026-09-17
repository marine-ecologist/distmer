# Derive estimable hypotheses for an unbalanced fixed model

Constructs sequential Type I, marginal Type II, or full-model Type III
projection hypotheses. Type II is restricted to terms without containing
higher-order interactions; Type III requires a full-rank model and
explicit sum-to-zero contrasts when interactions are present.

## Usage

``` r
derive_unbalanced_hypotheses(
  design,
  hypothesis = c("I", "II", "III"),
  terms = NULL,
  tol = 1e-10
)
```

## Arguments

- design:

  Canonical fixed-effect Design Object.

- hypothesis:

  One of `"I"`, `"II"`, or `"III"`.

- terms:

  Optional fixed terms to derive.

- tol:

  Numerical tolerance.

## Value

A named list of hypothesis projectors and estimability metadata.
