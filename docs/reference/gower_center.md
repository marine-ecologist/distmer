# Gower-centre a distance matrix

Implements `G = -0.5 * J %*% D^2 %*% J`, where `J = I - 11'/N`.

## Usage

``` r
gower_center(distance)
```

## Arguments

- distance:

  A `dist` object or full distance matrix.

## Value

The centred inner-product matrix `G`.
