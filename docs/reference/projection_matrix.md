# Orthogonal projection matrix

Uses a rank-revealing QR decomposition and therefore does not require an
explicit matrix inverse.

## Usage

``` r
projection_matrix(X, tol = 1e-10)
```

## Arguments

- X:

  Numeric design matrix.

- tol:

  QR tolerance.

## Value

A list containing the projection matrix and estimated rank.
