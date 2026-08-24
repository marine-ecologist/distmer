# 5. Pairwise contrasts and multiple comparisons

## comparison methods

When you compare specific groups — say treatment vs control — there are
two reasonable questions, and
[`permanova_contrasts()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_contrasts.md)
lets you choose which one you’re asking with the method argument.

`method = "isolated"` (the default) asks: “Do these two groups differ,
considering just those two groups?” It pulls out only the samples from
the levels being compared and runs a fresh two-group PERMANOVA on them —
its own error term, its own permutations (the approach that `PRIMER`
takes for pairwise comparisons). It makes no assumptions about the other
groups in your design.

`method = "pooled"` asks: “Do these two groups differ, using the error
estimated from the whole experiment?” It tests the same difference but
borrows the error term from the full model, which uses all your data and
so can have more power. The trade-off is that it assumes the groups are
similarly variable (homogeneous dispersion) — if one group is much more
scattered than another, this assumption is questionable.

Use `isolated` when you want a clean, assumption-light comparison that
matches `PRIMER` (i.e. default). Use `pooled` when your groups are
similarly variable and you want the extra power from the shared error
term.

## pairwise and multicomparisons:

### i) one-way pairwise contrasts:

Using the `dune` data: which management regimes differ in dune plant
composition?

all pairwise comparisons among the four management levels (6 contrasts),
Holm-adjusted:

``` r


data(dune, dune.env, package = "vegan")
bray <- vegdist(dune, method = "bray")
 
fit1 <- permanova(bray ~ Management, data = dune.env,
                  permutations = 9999, seed = 101)
summary(fit1)
#> Distance-based PERMANOVA
#> Engine: fixed effects  (Type I SS)
#> 
#> Factors
#>  Name       Type  Levels
#>  Management Fixed 4     
#> 
#> PERMANOVA table of results
#>  Source     df SS    MS    Pseudo-F P(perm) P(MC)  Unique perms
#>  Management 3  1.469 0.490 2.7672   0.0027  0.0054 9999        
#>  Residual   16 2.830 0.177                                     
#>  Total      19 4.299                                           
#> 
#> Details of the expected mean squares (EMS) for the model
#>  Source     EMS                           
#>  Management 1*V(Res) + 4.900*S(Management)
#>  Res        1*V(Res)                      
#> 
#> Construction of Pseudo-F ratio(s) from mean squares
#>  Source     Numerator    Denominator Num.df Den.df
#>  Management 1*Management 1*Res       3      16    
#> 
#> Estimates of components of variation
#>  Source        Estimate Sq.root
#>  S(Management) 0.064    0.253  
#>  V(Res)        0.177    0.421
 
permanova_contrasts(fit1, type = "pairwise") |> tab()
```

| contrast | pseudo_F | df_num | df_den | denominator | p_perm | unique_perms | estimable | p_adj |
|:---|---:|---:|---:|:---|---:|---:|:---|---:|
| BF - HF | 1.568 | 1 | 6 | MS\[full-model Residual\] | 0.200 | 9999 | TRUE | 0.241 |
| BF - NM | 2.715 | 1 | 7 | MS\[full-model Residual\] | 0.022 | 9999 | TRUE | 0.112 |
| BF - SF | 2.515 | 1 | 7 | MS\[full-model Residual\] | 0.059 | 9999 | TRUE | 0.178 |
| HF - NM | 3.423 | 1 | 9 | MS\[full-model Residual\] | 0.029 | 9999 | TRUE | 0.116 |
| HF - SF | 1.857 | 1 | 9 | MS\[full-model Residual\] | 0.121 | 9999 | TRUE | 0.241 |
| NM - SF | 3.426 | 1 | 10 | MS\[full-model Residual\] | 0.005 | 9999 | TRUE | 0.031 |

same, but control the family-wise error differently
(Benjamini-Hochberg):

``` r

permanova_contrasts(fit1, type = "pairwise", adjust = "BH")  |> tab()
```

| contrast | pseudo_F | df_num | df_den | denominator | p_perm | unique_perms | estimable | p_adj |
|:---|---:|---:|---:|:---|---:|---:|:---|---:|
| BF - HF | 1.568 | 1 | 6 | MS\[full-model Residual\] | 0.200 | 9999 | TRUE | 0.200 |
| BF - NM | 2.715 | 1 | 7 | MS\[full-model Residual\] | 0.022 | 9999 | TRUE | 0.058 |
| BF - SF | 2.515 | 1 | 7 | MS\[full-model Residual\] | 0.059 | 9999 | TRUE | 0.089 |
| HF - NM | 3.423 | 1 | 9 | MS\[full-model Residual\] | 0.029 | 9999 | TRUE | 0.058 |
| HF - SF | 1.857 | 1 | 9 | MS\[full-model Residual\] | 0.121 | 9999 | TRUE | 0.145 |
| NM - SF | 3.426 | 1 | 10 | MS\[full-model Residual\] | 0.005 | 9999 | TRUE | 0.031 |

### ii) treatment vs control:

Using the `pyrifos` dataset (12 ditches, 11 weeks, dose applied per
ditch): does each pesticide dose differ from the untreated control?

every dose vs control (dose = 0): 4 contrasts:

``` r

data(pyrifos, package = "vegan")
ditch <- gl(12, 1, length = 132)
week  <- gl(11, 12, labels = c(-4, -1, 0.1, 1, 2, 4, 8, 12, 15, 19, 24))
dose  <- factor(rep(c(0.1, 0, 0, 0.9, 0, 44, 6, 0.1, 44, 0.9, 0, 6), 11))
pyr.env  <- data.frame(ditch, week, dose)
pyr.bray <- vegdist(pyrifos)
 
fit2 <- permanova(pyr.bray ~ dose, data = pyr.env,
                  permutations = 999, seed = 101)
 

permanova_contrasts(fit2, type = "trt.vs.ctrl", ref = "0") |> tab()
```

| contrast | pseudo_F | df_num | df_den | denominator | p_perm | unique_perms | estimable | p_adj |
|:---|---:|---:|---:|:---|---:|---:|:---|---:|
| 0.1 - 0 | 2.313 | 1 | 64 | MS\[full-model Residual\] | 0.007 | 999 | TRUE | 0.007 |
| 0.9 - 0 | 3.963 | 1 | 64 | MS\[full-model Residual\] | 0.002 | 999 | TRUE | 0.004 |
| 6 - 0 | 6.224 | 1 | 64 | MS\[full-model Residual\] | 0.001 | 999 | TRUE | 0.004 |
| 44 - 0 | 10.726 | 1 | 64 | MS\[full-model Residual\] | 0.001 | 999 | TRUE | 0.004 |

### iii) polynomial contrasts

is there a monotone (linear) compositional trend across dose, and any
curvature?

``` r


permanova_contrasts(fit2, type = "poly", scores = c(0, 0.1, 0.9, 6, 44))  |> tab()
```

| contrast | pseudo_F | df_num | df_den | denominator | p_perm | unique_perms | estimable | p_adj |
|:---|---:|---:|---:|:---|---:|---:|:---|---:|
| linear | 10.781 | 1 | 127 | Residual (pooled; isolated undefined for \>2-level contrast) | NA | NA | TRUE | NA |
| quadratic | 5.018 | 1 | 127 | Residual (pooled; isolated undefined for \>2-level contrast) | NA | NA | TRUE | NA |
| cubic | 3.666 | 1 | 127 | Residual (pooled; isolated undefined for \>2-level contrast) | NA | NA | TRUE | NA |
| quartic | 2.854 | 1 | 127 | Residual (pooled; isolated undefined for \>2-level contrast) | NA | NA | TRUE | NA |

if `dose` were an ordered factor with numeric labels, scores are
inferred:

``` r

pyr.env$dose_ord <- ordered(pyr.env$dose,
                            levels = c("0", "0.1", "0.9", "6", "44"))
fit2b <- permanova(pyr.bray ~ dose_ord, data = pyr.env,
                   permutations = 999, seed = 101)
permanova_contrasts(fit2b, type = "poly")    |> tab()
```

| contrast | pseudo_F | df_num | df_den | denominator | p_perm | unique_perms | estimable | p_adj |
|:---|---:|---:|---:|:---|---:|---:|:---|---:|
| linear | 10.781 | 1 | 127 | Residual (pooled; isolated undefined for \>2-level contrast) | NA | NA | TRUE | NA |
| quadratic | 5.018 | 1 | 127 | Residual (pooled; isolated undefined for \>2-level contrast) | NA | NA | TRUE | NA |
| cubic | 3.666 | 1 | 127 | Residual (pooled; isolated undefined for \>2-level contrast) | NA | NA | TRUE | NA |
| quartic | 2.854 | 1 | 127 | Residual (pooled; isolated undefined for \>2-level contrast) | NA | NA | TRUE | NA |

### iii) sum to zero (grand mean)

``` r


permanova_contrasts(fit1, type = "contr.sum") |> tab()
```

| contrast | pseudo_F | df_num | df_den | denominator | p_perm | unique_perms | estimable | p_adj |
|:---|---:|---:|---:|:---|---:|---:|:---|---:|
| BF - mean | 2.515 | 1 | 7 | MS\[full-model Residual\] | 0.059 | 9999 | TRUE | 0.119 |
| HF - mean | 1.857 | 1 | 9 | MS\[full-model Residual\] | 0.121 | 9999 | TRUE | 0.121 |
| NM - mean | 3.426 | 1 | 10 | MS\[full-model Residual\] | 0.005 | 9999 | TRUE | 0.016 |

### iv) helmert contrasts

helmert (i.e. each level vs the mean of subsequent levels):

``` r


permanova_contrasts(fit1, type = "helmert") |> tab()
```

| contrast | pseudo_F | df_num | df_den | denominator | p_perm | unique_perms | estimable | p_adj |
|:---|---:|---:|---:|:---|---:|---:|:---|---:|
| helmert_1 | 1.568 | 1 | 6 | MS\[full-model Residual\] | 0.2 | 9999 | TRUE | 0.2 |
| helmert_2 | 4.761 | 1 | 16 | Residual (pooled; isolated undefined for \>2-level contrast) | NA | NA | TRUE | NA |
| helmert_3 | 2.829 | 1 | 16 | Residual (pooled; isolated undefined for \>2-level contrast) | NA | NA | TRUE | NA |

### v) custom contrasts

Test a specific ecological hypothesis: do the low doses (0, 0.1, 0.9)
together differ from the high doses (6, 44)? Columns MUST be in the
factor’s level order; check with levels():

``` r


Lm <- rbind(
  "low vs high"   = c( 1/3,  1/3,  1/3, -1/2, -1/2),    # note level order!
  "44 vs 6"       = c(   0,    0,    0,    1,   -1)
)
permanova_contrasts(fit2, type = "custom", custom = Lm) |> tab()
```

| contrast | pseudo_F | df_num | df_den | denominator | p_perm | unique_perms | estimable | p_adj |
|:---|---:|---:|---:|:---|---:|---:|:---|---:|
| low vs high | 11.754 | 1 | 127 | Residual (pooled; isolated undefined for \>2-level contrast) | NA | NA | TRUE | NA |
| 44 vs 6 | 2.536 | 1 | 42 | MS\[full-model Residual\] | 0.008 | 999 | TRUE | 0.008 |

a named list works too, one vector per contrast:

``` r


Llist <- list(
  "any dose vs control" = c(-4, 1, 1, 1, 1),            # control vs the rest
  "top dose vs control" = c(-1, 0, 0, 1, 0)             # 44 vs 0
)
permanova_contrasts(fit2, type = "custom", custom = Llist) |> tab()
```

| contrast | pseudo_F | df_num | df_den | denominator | p_perm | unique_perms | estimable | p_adj |
|:---|---:|---:|---:|:---|---:|---:|:---|---:|
| any dose vs control | 7.097 | 1 | 127 | Residual (pooled; isolated undefined for \>2-level contrast) | NA | NA | TRUE | NA |
| top dose vs control | 6.224 | 1 | 64 | MS\[full-model Residual\] | 0.001 | 999 | TRUE | 0.001 |

### vi) contrasts with random effects:

does treatment differ, tested at the correct level, with repeated
samples per site?

``` r


set.seed(7)
site  <- factor(rep(paste0("S", 1:8), each = 4))
treat <- factor(rep(c("ctrl", "trt"), each = 16))   # S1-S4 ctrl, S5-S8 trt
comm  <- matrix(rpois(32 * 15, 5), 32, 15)
d     <- vegdist(comm)
env   <- data.frame(site, treat)

fit3 <- permanova(d ~ treat + (1 | site), data = env,
                  permutations = 999, seed = 202)
 

permanova_contrasts(fit3, type = "pairwise") |> tab()
```

| contrast | pseudo_F | df_num | df_den | denominator | p_perm | unique_perms | estimable | p_adj |
|:---|---:|---:|---:|:---|---:|---:|:---|---:|
| ctrl - trt | 1.751 | 1 | 6 | MS\[site within fixed model\] | 0.115 | 999 | TRUE | 0.115 |

### contrast matrix

[`contrast_matrix()`](https://marine-ecologist.github.io/distmer/dev/reference/contrast_matrix.md)
is exported so you can see/verify the L matrix a family produces —
useful for checking custom coding or level order before you run
thousands of permutations.

``` r

contrast_matrix("trt.vs.ctrl", levels = levels(pyr.env$dose), ref = "0") 
#>          0 0.1 0.9 6 44
#> 0.1 - 0 -1   1   0 0  0
#> 0.9 - 0 -1   0   1 0  0
#> 6 - 0   -1   0   0 1  0
#> 44 - 0  -1   0   0 0  1
# contrast_matrix("poly", levels = levels(pyr.env$dose),
#                 scores = c(0, 0.1, 0.9, 44, 6))          # matches level order
# contrast_matrix("pairwise", levels = levels(dune.env$Management))
```

## programmatic use

``` r


library(distmer)
library(vegan)

data(pyrifos)

ditch <- gl(12, 1, length = 132)
week  <- gl(11, 12, labels = c(-4, -1, 0.1, 1, 2, 4, 8, 12, 15, 19, 24))
dose  <- factor(rep(c(0.1, 0, 0, 0.9, 0, 44, 6, 0.1, 44, 0.9, 0, 6), 11))

pyr.env <- data.frame(ditch, week, dose)
pyr.bray <- vegdist(pyrifos)  

fit <- permanova(
  pyr.bray ~ dose,
  data = pyr.env,
  permutations = 999,
  seed = 101
)

res <- permanova_contrasts(fit, type = "trt.vs.ctrl", ref = "0")
as.data.frame(res) |> tab()                      # full numeric table
summary(res)                             # same data.frame
subset(as.data.frame(res), p_adj < 0.05) # significant contrasts only
res$L         
```

\`\`\`
