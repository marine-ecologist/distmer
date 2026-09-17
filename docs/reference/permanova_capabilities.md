# Report the validated inference support matrix

Returns one row for every primary inferential engine in the validated
support contract. The registry is descriptive: it does not relax the
runtime design validators, which remain authoritative for a particular
compiled design.

## Usage

``` r
permanova_capabilities()
```

## Value

A data frame containing design envelopes, inferential scales,
denominators, permutation contracts, error-control claims, and
exclusions.
