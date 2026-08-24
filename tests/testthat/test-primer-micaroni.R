# Nested strata whose expected mean squares coincide. No data needed: PRIMER
# prints the expected mean squares and the denominator it builds from them.

source(testthat::test_path("fixtures", "primer-micaroni.R"))

test_that("a nested stratum matching the residual is still the denominator", {
  strata <- primer_micaroni_strata
  candidate <- do.call(rbind, primer_micaroni_ems[strata])
  rownames(candidate) <- strata
  # Mean squares consistent with the printed pseudo-F ratios; only their ratios
  # and the degrees of freedom matter to the construction.
  ms <- c(`ch(tr)` = 2.0, `sp(ch(tr))` = 1.3775, `ch(tr)xt` = 0.9, Res = 0.6)
  df <- c(`ch(tr)` = 4, `sp(ch(tr))` = 24, `ch(tr)xt` = 48, Res = 288)

  for (i in seq_len(nrow(primer_micaroni_construction))) {
    source_name <- primer_micaroni_construction$source[[i]]
    numerator <- primer_micaroni_ems[[source_name]]
    # Testing a random term sets its own component to zero.
    own <- paste0("V(", source_name, ")")
    if (own %in% names(numerator)) numerator[[own]] <- 0

    built <- .ems_f_construction(
      numerator, numerator_ms = 1, numerator_df =
        primer_micaroni_construction$num_df[[i]],
      candidate, ms, df, residual = "V(Res)", after = source_name
    )
    expect_false(is.null(built), info = source_name)
    expect_equal(built$denominator_df,
                 primer_micaroni_construction$den_df[[i]],
                 info = source_name)
    expect_equal(built$strata_used,
                 primer_micaroni_construction$denominator[[i]],
                 info = source_name)
  }
})

test_that("without the hierarchy the tie falls to the residual", {
  strata <- primer_micaroni_strata
  candidate <- do.call(rbind, primer_micaroni_ems[strata])
  rownames(candidate) <- strata
  ms <- c(`ch(tr)` = 2.0, `sp(ch(tr))` = 1.3775, `ch(tr)xt` = 0.9, Res = 0.6)
  df <- c(`ch(tr)` = 4, `sp(ch(tr))` = 24, `ch(tr)xt` = 48, Res = 288)
  numerator <- primer_micaroni_ems[["ch(tr)"]]
  numerator[["V(ch(tr))"]] <- 0

  # `after` is what carries the hierarchy. Without it the expected mean squares
  # are genuinely ambiguous -- sponge and residual read the same -- and the
  # residual is the only defensible choice.
  plain <- .ems_f_construction(numerator, 1, 4, candidate, ms, df,
                               residual = "V(Res)")
  expect_equal(plain$denominator_df, 288)
  guided <- .ems_f_construction(numerator, 1, 4, candidate, ms, df,
                                residual = "V(Res)", after = "ch(tr)")
  expect_equal(guided$denominator_df, 24)
})
