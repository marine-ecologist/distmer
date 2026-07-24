# Diagnose the Euclidean spectrum of a dissimilarity

Negative eigenvalues are retained and summarized rather than silently
discarded. Corrections, when added, will be separate sensitivity
analyses.

## Usage

``` r
distance_diagnostics(distance, tol = 1e-10)
```

## Arguments

- distance:

  A \`dist\` object or full distance matrix.

- tol:

  Eigenvalue tolerance.

## Value

A list containing the centred matrix, complete eigenvalue spectrum,
inertia summaries, and a Euclidean flag.
