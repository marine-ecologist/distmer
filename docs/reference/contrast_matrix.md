# Build a contrast (L) matrix from a named family

Build a contrast (L) matrix from a named family

## Usage

``` r
contrast_matrix(type, levels, ref = NULL, scores = NULL, custom = NULL)
```

## Arguments

- type:

  One of "pairwise","trt.vs.ctrl","poly","contr.sum","helmert",
  "custom", or a matrix / named list of contrast rows.

- levels:

  Character vector of factor levels, in model order.

- ref:

  Reference level for "trt.vs.ctrl".

- scores:

  Optional numeric spacing for "poly".

- custom:

  Matrix or named list of rows when type = "custom".

## Value

Numeric matrix: rows = contrasts (named), columns = levels (named).
