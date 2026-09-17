# Derive fixed-model marginal means

Constructs coefficient-space rows for every level of a target factor in
a full-rank categorical fixed model. Equal weighting averages the
complete nuisance-factor reference grid. Proportional weighting averages
observed nuisance combinations by their joint frequency.

## Usage

``` r
derive_marginal_means(
  design,
  factor,
  weights = c("equal", "proportional"),
  tol = 1e-10
)
```

## Arguments

- design:

  Canonical fixed-effect Design Object.

- factor:

  Target categorical predictor.

- weights:

  `"equal"` or `"proportional"` nuisance-cell weights.

- tol:

  Numerical tolerance.

## Value

Marginal-mean coefficient rows and reference-grid metadata.
