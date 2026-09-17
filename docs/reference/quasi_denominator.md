# Construct an EMS-matched quasi-denominator

Finds coefficients for candidate mean squares whose combined EMS matches
the numerator EMS. An exact single mean square is preferred; otherwise a
minimum-norm linear combination is returned with Satterthwaite degrees
of freedom.

## Usage

``` r
quasi_denominator(
  numerator_ems,
  denominator_ems,
  mean_squares,
  dfs,
  tol = 1e-10
)
```

## Arguments

- numerator_ems:

  Named numeric EMS coefficient vector.

- denominator_ems:

  Matrix with one candidate mean square per row and EMS components in
  columns.

- mean_squares:

  Candidate observed mean squares.

- dfs:

  Candidate degrees of freedom.

- tol:

  Numerical matching tolerance.

## Value

A machine-readable denominator specification.
