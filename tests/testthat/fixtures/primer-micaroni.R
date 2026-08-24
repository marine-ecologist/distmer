# PRIMER 7 / PERMANOVA+ output from Micaroni et al., sponge deoxygenation
# experiment: treatment fixed, chamber random within treatment, sponge random
# within chamber, time fixed, Type I sums of squares, restricted ("fixed effects
# sum to zero") expected mean squares.
#
# The point of this model for distmer is the sponge stratum. PRIMER prints its
# expected mean square as 1*V(Res) -- identical to the residual's -- so the
# expected mean squares alone cannot say which of the two is the denominator for
# chamber. The nesting hierarchy decides it: chamber is tested over sponge, not
# over the residual, giving 24 denominator df rather than 288.

primer_micaroni_ems <- list(
  tr           = c(`V(Res)` = 1, `V(ch(tr))` = 65, `V(ch(tr)xt)` = 0),
  `ch(tr)`     = c(`V(Res)` = 1, `V(ch(tr))` = 65, `V(ch(tr)xt)` = 0),
  `sp(ch(tr))` = c(`V(Res)` = 1, `V(ch(tr))` = 0,  `V(ch(tr)xt)` = 0),
  t            = c(`V(Res)` = 1, `V(ch(tr))` = 0,  `V(ch(tr)xt)` = 5),
  trxt         = c(`V(Res)` = 1, `V(ch(tr))` = 0,  `V(ch(tr)xt)` = 5),
  `ch(tr)xt`   = c(`V(Res)` = 1, `V(ch(tr))` = 0,  `V(ch(tr)xt)` = 5),
  Res          = c(`V(Res)` = 1, `V(ch(tr))` = 0,  `V(ch(tr)xt)` = 0)
)

# Source, df and the denominator PRIMER builds, from the workspace.
primer_micaroni_construction <- data.frame(
  source      = c("tr", "ch(tr)", "sp(ch(tr))", "t", "trxt", "ch(tr)xt"),
  denominator = c("ch(tr)", "sp(ch(tr))", "Res", "ch(tr)xt", "ch(tr)xt", "Res"),
  num_df      = c(1, 4, 24, 12, 12, 48),
  den_df      = c(4, 24, 288, 48, 48, 288),
  stringsAsFactors = FALSE
)

# Strata run coarsest to finest, as PRIMER prints them.
primer_micaroni_strata <- c("ch(tr)", "sp(ch(tr))", "ch(tr)xt", "Res")
