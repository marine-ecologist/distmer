# 5. Pairwise contrasts and multiple comparisons

An overall model can show that community composition differs among
groups, but it does not say which ecological comparisons are
responsible. This guide uses the four management classes in
[`vegan::dune`](https://vegandevs.github.io/vegan/reference/dune.html)
to move from one overall question to planned contrasts, all pairwise
comparisons, and multiplicity control.

``` r

management <- data.frame(Management = droplevels(dune.env$Management))
m1 <- permanova(bray ~ Management, management, permutations = 499, seed = 41)
summary(m1) |> 
  kbl(digits=3) |> 
  kable_styling(bootstrap_options = c("striped", "hover", "condensed")) 
#> Distance-based PERMANOVA
#> Engine: fixed effects  (Type I SS)
#> 
#> Factors
#>  Name       Type  Levels
#>  Management Fixed 4     
#> 
#> PERMANOVA table of results
#>  Source     df SS    MS    Pseudo-F P(perm) P(MC)  Unique perms
#>  Management 3  1.469 0.490 2.7672   0.0040  0.0050 499         
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
```

| Source     |  df |    SS |    MS | Pseudo-F | P(perm) | P(MC) | Unique perms |
|:-----------|----:|------:|------:|---------:|--------:|------:|-------------:|
| Management |   3 | 1.469 | 0.490 |    2.767 |   0.004 | 0.005 |          499 |
| Residual   |  16 | 2.830 | 0.177 |       NA |      NA |    NA |           NA |
| Total      |  19 | 4.299 |    NA |       NA |      NA |    NA |           NA |

The ordination is a visual summary, not the test itself. The model uses
the complete Bray–Curtis matrix, whereas the plot shows only two axes.

``` r

ord <- cmdscale(bray, k = 2)
ord_plot <- data.frame(
  Axis1 = ord[, 1], Axis2 = ord[, 2], Management = management$Management
)

ggplot() +
  theme_bw() + 
  geom_point(data = ord_plot, aes(Axis1, Axis2, fill = Management), shape=21, size = 3) +
  stat_ellipse(data = ord_plot, aes(Axis1, Axis2, colour = Management), linewidth = 0.6, type = "norm",
               show.legend = FALSE) +
  labs(title = "Dune communities by management",
       x = "PCoA axis 1", y = "PCoA axis 2", colour = "Management")
#> Too few points to calculate an ellipse
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_path()`).
```

![Dune sites in a two-axis PCoA coloured by management
class](contrasts_files/figure-html/ordination-1.png)

## Start with ecological questions

A planned contrast is most useful when it expresses a question chosen
before looking at the pairwise results. Here we compare biological
farming with hobby farming, and the average of the two more intensively
managed classes with the average of the other two classes. Coefficients
are applied to equally weighted group centroids and must sum to zero.

``` r

design <- compile_permanova_design(~ Management, management, bray)
planned <- rbind(
  `biological - hobby farming` = c(BF = 1, HF = -1, NM = 0, SF = 0),
  `BF/HF mean - NM/SF mean` = c(BF = 0.5, HF = 0.5, NM = -0.5, SF = -0.5)
)

planned_fit <- permanova_contrasts(
  design, planned, adjust = "holm", permutations = 499, seed = 41
)
planned_fit$table |> 
  kbl(digits=3) |> 
  kable_styling(bootstrap_options = c("striped", "hover", "condensed")) 
```

|  | contrast | df | ss | ms | denominator_df | pseudo_f | p_perm | p_adjusted | adjustment |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|:---|
| biological - hobby farming | biological - hobby farming | 1 | 0.162 | 0.162 | 16 | 0.914 | 0.458 | 0.458 | holm |
| BF/HF mean - NM/SF mean | BF/HF mean - NM/SF mean | 1 | 0.563 | 0.563 | 16 | 3.181 | 0.020 | 0.040 | holm |

Each contrast keeps the remaining group structure in the model, then
obtains a permutation p-value for that comparison. Holm adjustment
protects the family of planned questions without assuming that the
contrasts are independent.

## All pairwise comparisons

When there is no small set of pre-specified comparisons,
[`permanova_pairwise()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_pairwise.md)
constructs every pair. The raw p-value describes one comparison; the
adjusted value accounts for having examined the complete family.

``` r

pw <- permanova_pairwise(
  design, adjust = "holm", permutations = 499, seed = 41
)
pw$table[, c("contrast", "pseudo_f", "p_perm", "p_adjusted")] |> 
  kbl(digits=3) |> 
  kable_styling(bootstrap_options = c("striped", "hover", "condensed")) 
```

|         | contrast | pseudo_f | p_perm | p_adjusted |
|:--------|:---------|---------:|-------:|-----------:|
| BF - HF | BF - HF  |    0.914 |  0.458 |      0.458 |
| BF - NM | BF - NM  |    3.201 |  0.012 |      0.048 |
| BF - SF | BF - SF  |    2.271 |  0.076 |      0.228 |
| HF - NM | HF - NM  |    3.682 |  0.006 |      0.030 |
| HF - SF | HF - SF  |    1.599 |  0.150 |      0.300 |
| NM - SF | NM - SF  |    4.282 |  0.002 |      0.012 |

``` r

pw_plot <- pw$table
pw_plot$contrast <- reorder(pw_plot$contrast, pw_plot$p_adjusted)

ggplot(pw_plot, aes(p_adjusted, contrast)) +
  geom_vline(xintercept = 0.05, linetype = 2, colour = "grey45") +
  geom_segment(aes(x = 0, xend = p_adjusted, yend = contrast),
               colour = "grey75") +
  geom_point(size = 2.7, colour = "#0072B2") +
  scale_x_continuous(limits = c(0, 1)) +
  labs(title = "Pairwise management comparisons",
       x = "Holm-adjusted permutation p-value", y = NULL)
```

![Adjusted p-values for all pairwise management
comparisons](contrasts_files/figure-html/pairwise-plot-1.png)

Do not interpret a non-significant pair as proof that its communities
are equivalent. It says only that this analysis did not resolve a
difference at the chosen error rate and permutation depth.

## Comparisons adjusted for other predictors

In a multifactor model, a management contrast needs an explicit rule for
averaging over other categorical predictors. `weights = "equal"`
describes an equally represented ecological population;
`weights = "proportional"` describes the observed sample composition.

``` r

factors <- data.frame(
  Management = droplevels(dune.env$Management),
  Use = droplevels(dune.env$Use)
)
multifactor <- compile_permanova_design(~ Management + Use, factors, bray)

permanova_marginal_contrasts(
  multifactor, factor = "Management", contrasts = planned,
  weights = "equal", permutations = 199, seed = 41
)$table |> 
  kbl(digits=3) |> 
  kable_styling(bootstrap_options = c("striped", "hover", "condensed")) 
```

|  | contrast | weighting | df | ss | denominator_df | pseudo_f | p_perm | p_adjusted | adjustment |
|:---|:---|:---|---:|---:|---:|---:|---:|---:|:---|
| biological - hobby farming | biological - hobby farming | equal | 1 | 0.159 | 14 | 0.908 | 0.495 | 0.495 | holm |
| BF/HF mean - NM/SF mean | BF/HF mean - NM/SF mean | equal | 1 | 0.652 | 14 | 3.721 | 0.015 | 0.030 | holm |

If a factor interacts with a numeric environmental gradient, its
comparison also depends on where that gradient is evaluated. The
reference value is therefore supplied rather than silently chosen.

``` r

env <- data.frame(
  Management = droplevels(dune.env$Management),
  Use = droplevels(dune.env$Use),
  A1 = dune.env$A1
)
adjusted_design <- compile_permanova_design(
  ~ Management * A1 + Use, env, bray
)

permanova_adjusted_contrasts(
  adjusted_design, factor = "Management", contrasts = planned,
  at = c(A1 = mean(env$A1)), weights = "equal",
  permutations = 199, seed = 41
)$table |> 
  kbl(digits=3) |> 
  kable_styling(bootstrap_options = c("striped", "hover", "condensed")) 
```

|  | contrast | weighting | df | ss | denominator_df | pseudo_f | p_perm | p_adjusted | adjustment |
|:---|:---|:---|---:|---:|---:|---:|---:|---:|:---|
| biological - hobby farming | biological - hobby farming | equal | 1 | 0.100 | 10 | 0.628 | 0.64 | 1 | holm |
| BF/HF mean - NM/SF mean | BF/HF mean - NM/SF mean | equal | 1 | 0.095 | 10 | 0.599 | 0.67 | 1 | holm |

## Choosing familywise protection

For most ecological reports, Holm-adjusted planned or pairwise
comparisons are a clear default. `distmer` also provides two
shared-permutation procedures:

- [`permanova_maxT_pairwise()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_maxT_pairwise.md)
  compares every statistic with the largest value in each common-null
  permutation. Its familywise guarantee is for the complete null in
  which all group centroids coincide.
- [`permanova_pairwise_strong()`](https://marine-ecologist.github.io/distmer/dev/reference/permanova_pairwise_strong.md)
  uses closed testing and protects every pair even when other pairs
  genuinely differ.

``` r

max_f <- permanova_maxT_pairwise(design, permutations = 199, seed = 41)
closed <- permanova_pairwise_strong(design, permutations = 199, seed = 41)

max_f$table[, c("contrast", "p_perm", "p_maxT")] |> 
  kbl(digits=3) |> 
  kable_styling(bootstrap_options = c("striped", "hover", "condensed")) 
```

|         | contrast | p_perm | p_maxT |
|:--------|:---------|-------:|-------:|
| BF - HF | BF - HF  |  0.495 |  0.910 |
| BF - NM | BF - NM  |  0.005 |  0.065 |
| BF - SF | BF - SF  |  0.060 |  0.190 |
| HF - NM | HF - NM  |  0.010 |  0.025 |
| HF - SF | HF - SF  |  0.115 |  0.465 |
| NM - SF | NM - SF  |  0.010 |  0.010 |

``` r

closed$table[, c("contrast", "p_perm", "p_closed")] |> 
  kbl(digits=3) |> 
  kable_styling(bootstrap_options = c("striped", "hover", "condensed")) 
```

|         | contrast | p_perm | p_closed |
|:--------|:---------|-------:|---------:|
| BF - HF | BF - HF  |  0.196 |    0.196 |
| BF - NM | BF - NM  |  0.024 |    0.075 |
| BF - SF | BF - SF  |  0.060 |    0.075 |
| HF - NM | HF - NM  |  0.035 |    0.060 |
| HF - SF | HF - SF  |  0.095 |    0.170 |
| NM - SF | NM - SF  |  0.015 |    0.025 |

The strong finite-sample guarantee has an important boundary. It applies
to pairwise equality of **complete distributions** for independent
one-factor observations. It is not a strong guarantee for centroid-only
nulls when groups have different dispersion or shape. Use the distance
and dispersion checks in
[`vignette("dispersion")`](https://marine-ecologist.github.io/distmer/dev/articles/dispersion.md)
alongside these comparisons.

## What to report

Report the contrast coefficients or pair labels, their ecological
rationale, the raw and adjusted p-values, the adjustment method, the
dissimilarity, the number of permutations, and the seed. A pseudo-F is
an effect-to-residual ratio, not an effect size; show an ordination or
centroid-distance summary as well.
