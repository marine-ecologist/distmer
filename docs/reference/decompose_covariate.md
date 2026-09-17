# Split a covariate into between-unit and within-unit components

Separates one or more continuous covariates into a between-unit part
(the grouping-unit mean) and a within-unit part (the deviation from that
mean). Entering both parts as predictors lets engines that support both
strata adjust the between-unit and within-unit relationships separately.
The within part is a genuine observation-level covariate suitable for
[`permanova_nested_covariate()`](https://marine-ecologist.github.io/distmer/reference/permanova_nested_covariate.md);
the between part is a unit-level summary and must only be entered in an
engine that assigns numeric terms to that level (for example
[`permanova_mixed_covariate()`](https://marine-ecologist.github.io/distmer/reference/permanova_mixed_covariate.md)).

## Usage

``` r
decompose_covariate(
  data,
  covariate,
  by,
  center = TRUE,
  suffix = c("_between", "_within")
)
```

## Arguments

- data:

  A data frame.

- covariate:

  Character vector of continuous covariate column names.

- by:

  Character vector of grouping column name(s) defining the units; the
  between/within split is taken with respect to their combination.

- center:

  Grand-centre the between-unit component (default `TRUE`) so it is a
  deviation from the overall mean. The within-unit component always has
  mean zero within every unit regardless of this setting.

- suffix:

  Length-two character vector of column-name suffixes for the between
  and within parts.

## Value

`data` with two added columns per covariate. With `center = FALSE` the
between and within columns sum to the original covariate.

## Examples

``` r
df <- data.frame(site = factor(rep(1:3, each = 4)), x = rnorm(12))
decompose_covariate(df, "x", by = "site")
#>    site            x  x_between   x_within
#> 1     1 -1.400043517 -0.5353601 -0.5031532
#> 2     1  0.255317055 -0.5353601  1.1522074
#> 3     1 -2.437263611 -0.5353601 -1.5403733
#> 4     1 -0.005571287 -0.5353601  0.8913191
#> 5     2  0.621552721  0.2867355  0.6963474
#> 6     2  1.148411606  0.2867355  1.2232063
#> 7     2 -1.821817661  0.2867355 -1.7470230
#> 8     2 -0.247325302  0.2867355 -0.1725306
#> 9     3 -0.244199607  0.2486246 -0.1312940
#> 10    3 -0.282705449  0.2486246 -0.1697998
#> 11    3 -0.553699384  0.2486246 -0.4407938
#> 12    3  0.628982042  0.2486246  0.7418876
```
