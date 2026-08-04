# Run the deterministic validation suite

Run the deterministic validation suite

## Usage

``` r
run_permanova_validation(seed = 1L, include_vegan = TRUE)
```

## Arguments

- seed:

  Reproducibility seed.

- include_vegan:

  Compare the supported fixed one-way case with
  [`vegan::adonis2()`](https://vegandevs.github.io/vegan/reference/adonis.html)
  when vegan is installed.

## Value

A data frame of validation checks with class `permanova_validation`.
