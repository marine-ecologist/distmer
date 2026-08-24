# Pairwise contrasts within the cells of other factors

Contrasts the levels of one fixed factor separately inside every
combination of the `within` factors, as PRIMER's "within level" pairwise
tables do. Each pair gets its own denominator: the declared random
strata are rebuilt on that subset, any left with no degrees of freedom
are dropped, and the contrast's expected mean square is matched against
whichever survive. Three pairs inside one cell can therefore carry three
different denominators and three different Satterthwaite degrees of
freedom.

## Usage

``` r
permanova_contrasts_within(
  fit,
  factor,
  within,
  permutations = 999L,
  seed = NULL,
  adjust = "none",
  max_exact = 100000L,
  tol = 1e-10
)
```

## Arguments

- fit:

  A `permanova` object.

- factor:

  Name of the fixed factor whose levels are contrasted.

- within:

  Character vector of factors whose combinations define the cells.

- permutations:

  Permutations per contrast.

- seed:

  Optional random seed.

- adjust:

  Multiplicity method passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html), applied
  within each cell rather than across the whole table.

- max_exact:

  Maximum permutation space allowed for exact enumeration.

- tol:

  Numerical tolerance.

## Value

A data frame of contrasts, one row per pair per cell.
