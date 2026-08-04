# Structured (Cornfield–Tukey / Satterthwaite) denominator construction

Builds a pseudo-F denominator for a mixed-model term by *synthesising*
mean squares in containment-lattice order, rather than by the
unconstrained minimum-norm match used in
[`quasi_denominator()`](https://marine-ecologist.github.io/distmer/dev/reference/quasi_denominator.md).
Given the expected mean squares (EMS) of a set of sources, each random
source's variance component, and the observed mean squares, it peels the
target EMS one variance component at a time from the top of the lattice
down, using the source that directly estimates each component. On a
clean balanced design this returns the exact single-source denominator;
on a perturbed design (e.g. an unbalanced or covariate-containing mixed
model) it returns the structured linear combination with Satterthwaite
degrees of freedom, matching the construction used by classical
mixed-model software.

## Usage

``` r
structured_denominator(ems, mean_squares, dfs, leading, target, tol = 1e-08)
```

## Arguments

- ems:

  Numeric matrix of EMS coefficients, one row per source (named) and one
  column per variance component (named). Entry `ems[i, j]` is the
  coefficient of component `j` in the expected mean square of source
  `i`.

- mean_squares:

  Named numeric vector of observed mean squares, one per source (names
  matching `rownames(ems)`).

- dfs:

  Named numeric vector of degrees of freedom, one per source.

- leading:

  Named character vector mapping each *random* source to the variance
  component it directly estimates (e.g.
  `c(Residual = "residual", Site = "site")`); fixed sources are omitted
  or `NA`.

- target:

  The source to test: the denominator is built to match its EMS with its
  own leading component removed. May instead be a named numeric target
  EMS vector.

- tol:

  Numerical tolerance for treating coefficients as zero.

## Value

A list with the denominator `weights` (named), its value, and its
Satterthwaite `df`.

## Examples

``` r
# A balanced two-level nested design: Site over Residual is exact.
ems <- rbind(Site = c(residual = 1, site = 4),
             Residual = c(residual = 1, site = 0))
structured_denominator(ems, c(Site = 3, Residual = 1), c(Site = 5, Residual = 40),
                       leading = c(Residual = "residual", Site = "site"),
                       target = "Site")
#> $weights
#> Residual 
#>        1 
#> 
#> $value
#> [1] 1
#> 
#> $df
#> [1] 40
#> 
```
