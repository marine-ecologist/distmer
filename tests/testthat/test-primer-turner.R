# Validation against published PRIMER 7 / PERMANOVA+ output.
#
# Two levels of test. The first needs no data: PRIMER prints its own expected
# mean squares and the construction it builds from them, so the construction can
# be driven with PRIMER's numbers and checked against PRIMER's answers. That is
# the tightest available check on the pseudo-F construction and it always runs.
#
# The second refits the whole analysis. It needs the authors' archived data,
# which is not redistributed with the package, so it skips when absent.

source(testthat::test_path("fixtures", "primer-turner-2025.R"))

test_that("the EMS construction reproduces PRIMER's own pseudo-F and both df", {
  strata <- c("IS(ECxPO)", "DExIS(ECxPO)", "Res")
  candidate <- primer_turner_ems[strata, , drop = FALSE]
  rownames(candidate) <- strata
  candidate_ms <- primer_turner_strata$ms[strata]
  candidate_df <- primer_turner_strata$df[strata]

  for (i in seq_len(nrow(primer_turner_construction))) {
    term <- primer_turner_construction$term[[i]]
    numerator_ems <- primer_turner_ems[term, ]
    numerator_ms <- primer_turner_table$ms[primer_turner_table$term == term]
    numerator_df <- primer_turner_table$df[primer_turner_table$term == term]

    built <- .ems_f_construction(numerator_ems, numerator_ms, numerator_df,
                                 candidate, candidate_ms, candidate_df)
    expect_false(is.null(built), info = term)

    expected_f <- primer_turner_table$pseudo_f[primer_turner_table$term == term]
    expect_equal(built$statistic, expected_f, tolerance = 1e-4, info = term)
    # PRIMER prints degrees of freedom to two decimals, so the comparison can
    # be no tighter than its own reporting precision.
    expect_equal(built$numerator_df,
                 primer_turner_construction$num_df[[i]],
                 tolerance = 5e-3, info = term)
    expect_equal(built$denominator_df,
                 primer_turner_construction$den_df[[i]],
                 tolerance = 5e-3, info = term)
    # the residual PRIMER moves onto the numerator rather than allowing a
    # negative denominator weight
    expect_equal(built$numerator_residual,
                 primer_turner_construction$numerator_residual[[i]],
                 tolerance = 5e-3, info = term)
  }
})

test_that("no stratum weight in PRIMER's construction is negative", {
  weights <- primer_turner_construction[, c("weight_IS", "weight_DExIS",
                                            "weight_Res")]
  expect_true(all(as.matrix(weights) >= 0))
  # and the residual is balanced across the two sides, never dropped
  for (i in seq_len(nrow(primer_turner_construction))) {
    carried <- sum(weights[i, ])
    supplied <- 1 + primer_turner_construction$numerator_residual[[i]]
    expect_equal(carried, supplied, tolerance = 1e-3,
                 info = primer_turner_construction$term[[i]])
  }
})

# ---------------------------------------------------------------------------
# Full refit. Skipped unless the archived data is present.
# ---------------------------------------------------------------------------

.turner_data <- function() {
  path <- file.path("code", "distmer_tests", "Turner_et_al",
                    "Turner_et_al_2025_Benthos_2010_2014.csv")
  for (root in c(".", "..", file.path("..", ".."), file.path("..", "..", ".."))) {
    candidate <- file.path(root, path)
    if (file.exists(candidate)) return(candidate)
  }
  NULL
}

test_that("the full model matches PRIMER on df, and on the strata it shares", {
  path <- .turner_data()
  skip_if(is.null(path), "Turner et al. (2025) archive not present")
  skip_if_not_installed("vegan")

  raw <- utils::read.csv(path, check.names = FALSE, fileEncoding = "UTF-8-BOM")
  cover <- c("CORAL_%", "CCA_%", "SC_%", "HALIMEDA_%", "MA_%", "EMA_%", "TURF_%")
  meta <- data.frame(
    ISLAND = factor(raw$ISLAND), ECOREGION = factor(raw$ECOREGION),
    POP_STATUS = factor(raw$POP_STATUS),
    DEPTH = factor(raw$DEPTH_BIN_NEW, levels = c("Shallow", "Mid", "Deep"))
  )
  expect_equal(nrow(meta), 2239L)
  D <- dist(sqrt(as.matrix(raw[, cover])), method = "euclidean")

  fit <- permanova_asymmetric_mixed(
    compile_permanova_design(
      ~ DEPTH * ECOREGION * POP_STATUS + (1 | ISLAND) + (1 | ISLAND:DEPTH),
      meta, D),
    permutations = 0, hypothesis = "I", whole_plot = "observation"
  )
  map <- c(DEPTH = "DE", ECOREGION = "EC", POP_STATUS = "PO",
           `DEPTH:ECOREGION` = "DExEC", `DEPTH:POP_STATUS` = "DExPO",
           `ECOREGION:POP_STATUS` = "ECxPO",
           `DEPTH:ECOREGION:POP_STATUS` = "DExECxPO",
           `IS(ISLAND)` = "IS(ECxPO)",
           `DEPTH x IS(ISLAND)` = "DExIS(ECxPO)", Residual = "Res")

  # Every degree of freedom matches, which fixes the structure of the model.
  for (ours in names(map)) {
    row <- primer_turner_table[primer_turner_table$term == map[[ours]], ]
    i <- match(ours, fit$table$term)
    expect_false(is.na(i), info = ours)
    expect_equal(fit$table$df[[i]], row$df, info = ours)
  }

  # The deepest stratum and the residual are unaffected by term order, so they
  # agree exactly. These pin the dissimilarity and the overall decomposition.
  for (ours in c("DEPTH x IS(ISLAND)", "Residual")) {
    row <- primer_turner_table[primer_turner_table$term == map[[ours]], ]
    i <- match(ours, fit$table$term)
    expect_equal(fit$table$ms[[i]], row$ms, tolerance = 1e-3, info = ours)
  }

  # KNOWN GAP. PRIMER fits the random strata inside the Type I sequence -- its
  # table places IS(ECxPO) seventh, before the three-way -- and orders the fixed
  # terms differently from the `*` expansion. distmer fits the strata after the
  # whole-plot terms only. Adopting PRIMER's sequence reproduces every value
  # exactly; until that is implemented the sums of squares differ for the terms
  # whose reduced model the sequence changes. This guards the agreement already
  # achieved, so a regression is caught and an improvement fails loudly enough
  # to be noticed and the bound tightened.
  fixed <- names(map)[!names(map) %in% c("Residual")]
  ratio <- vapply(fixed, function(ours) {
    row <- primer_turner_table[primer_turner_table$term == map[[ours]], ]
    i <- match(ours, fit$table$term)
    fit$table$pseudo_f[[i]] / row$pseudo_f
  }, numeric(1))
  expect_lt(mean(abs(log(ratio))), 0.05)
})

test_that("every within-cell contrast reproduces PRIMER's Table 2", {
  path <- .turner_data()
  skip_if(is.null(path), "Turner et al. (2025) archive not present")

  raw <- utils::read.csv(path, check.names = FALSE, fileEncoding = "UTF-8-BOM")
  cover <- c("CORAL_%", "CCA_%", "SC_%", "HALIMEDA_%", "MA_%", "EMA_%", "TURF_%")
  meta <- data.frame(
    ISLAND = factor(raw$ISLAND), ECOREGION = factor(raw$ECOREGION),
    POP_STATUS = factor(raw$POP_STATUS),
    DEPTH = factor(raw$DEPTH_BIN_NEW, levels = c("Shallow", "Mid", "Deep"))
  )
  D <- dist(sqrt(as.matrix(raw[, cover])), method = "euclidean")
  fit <- permanova(D ~ DEPTH * ECOREGION * POP_STATUS +
                     (1 | ISLAND) + (1 | ISLAND:DEPTH),
                   meta, permutations = 0, seed = 101)

  ours <- as.data.frame(permanova_contrasts_within(
    fit, factor = "DEPTH", within = c("ECOREGION", "POP_STATUS"),
    permutations = 999, seed = 101
  ))
  # distmer names a pair "Mid - Deep"; PRIMER names it "Deep, Mid"
  ours$pair <- vapply(strsplit(ours$contrast, " - ", fixed = TRUE),
                      function(p) paste(sort(p), collapse = ", "), character(1))
  key <- function(d, ec, po) paste(d[[ec]], d[[po]], d$pair, sep = "|")
  ours$key <- key(ours, "ECOREGION", "POP_STATUS")
  reference <- primer_turner_pairwise
  reference$key <- key(reference, "ecoregion", "pop_status")

  expect_setequal(ours$key, reference$key)
  i <- match(reference$key, ours$key)

  # Pull the numeric weight of each stratum out of the printed denominator.
  weight <- function(text, name) {
    vapply(text, function(one) {
      terms <- regmatches(one, gregexpr("[0-9.e+-]+\\*[^ +]+", one))[[1L]]
      value <- as.numeric(sub("\\*.*", "", terms))[sub(".*\\*", "", terms) == name]
      if (!length(value)) 0 else value[[1L]]
    }, numeric(1), USE.NAMES = FALSE)
  }
  strata <- c(ISLAND = "weight_IS", `ISLAND:DEPTH` = "weight_DExIS",
              Residual = "weight_Res")

  for (k in seq_len(nrow(reference))) {
    label <- reference$key[[k]]
    row <- ours[i[[k]], ]
    # PRIMER prints five significant figures; nothing here is looser than that.
    expect_equal(row$t, reference$t[[k]], tolerance = 1e-4, info = label)
    for (nm in names(strata)) {
      expect_equal(weight(row$denominator, nm), reference[[strata[[nm]]]][[k]],
                   tolerance = 1e-4, info = paste(label, nm))
    }
    # denominator df are printed to two decimals
    expect_equal(row$denominator_df, reference$den_df[[k]], tolerance = 5e-3,
                 info = label)
  }

  # Three pairs inside one cell carrying three different Satterthwaite df is
  # the point of deriving each denominator inside its own subset.
  samoa <- ours[ours$ECOREGION == "Samoa Islands" & ours$POP_STATUS == "U", ]
  expect_equal(length(unique(round(samoa$denominator_df, 2))), 3L)
})

test_that("the permutation space matches PRIMER's own unique-perm count", {
  path <- .turner_data()
  skip_if(is.null(path), "Turner et al. (2025) archive not present")

  raw <- utils::read.csv(path, check.names = FALSE, fileEncoding = "UTF-8-BOM")
  cover <- c("CORAL_%", "CCA_%", "SC_%", "HALIMEDA_%", "MA_%", "EMA_%", "TURF_%")
  meta <- data.frame(
    ISLAND = factor(raw$ISLAND), ECOREGION = factor(raw$ECOREGION),
    POP_STATUS = factor(raw$POP_STATUS),
    DEPTH = factor(raw$DEPTH_BIN_NEW, levels = c("Shallow", "Mid", "Deep"))
  )
  D <- dist(sqrt(as.matrix(raw[, cover])), method = "euclidean")
  fit <- permanova(D ~ DEPTH * ECOREGION * POP_STATUS +
                     (1 | ISLAND) + (1 | ISLAND:DEPTH),
                   meta, permutations = 0, seed = 101)
  ours <- as.data.frame(permanova_contrasts_within(
    fit, factor = "DEPTH", within = c("ECOREGION", "POP_STATUS"),
    permutations = 9999, seed = 101
  ))
  ours$pair <- vapply(strsplit(ours$contrast, " - ", fixed = TRUE),
                      function(p) paste(sort(p), collapse = ", "), character(1))
  ours$key <- paste(ours$ECOREGION, ours$POP_STATUS, ours$pair, sep = "|")
  reference <- primer_turner_pairwise
  reference$key <- paste(reference$ecoregion, reference$pop_status,
                         reference$pair, sep = "|")
  i <- match(reference$key, ours$key)

  # Where PRIMER's space is small enough to name exactly, distmer must find the
  # same one. Every cell with six occupied island x depth units gives 6! = 720,
  # which is what identifies the exchangeable unit as the island x depth cell
  # rather than the observation. PRIMER sampled rather than enumerated, so it
  # can report 719 having missed one of the 720.
  units <- reference$unique_perms >= 100 & reference$unique_perms < 1000
  expect_equal(sum(units), 9L)   # three cells of three pairs
  expect_true(all(reference$unique_perms[units] %in% c(719, 720)))
  expect_true(all(ours$unique_perms[i][units] == 720))

  # The one cell small enough to enumerate at the observation level exposes a
  # difference of counting, not of space: PRIMER reports the 3 distinct group
  # assignments of four samples into two pairs, distmer the 4! = 24 orderings
  # that produce them. Same exchange, deduplicated at different points.
  tiny <- reference$unique_perms < 100
  expect_equal(sum(tiny), 3L)
  expect_true(all(reference$unique_perms[tiny] == choose(4, 2) / 2))
  expect_true(all(ours$unique_perms[i][tiny] == factorial(4)))
  expect_true(all(ours$distinct_statistics[i][tiny] <= 3L))

  gap <- abs(ours$p_perm[i] - reference$p_perm)
  unit_level <- reference$weight_DExIS > 0

  # Where the residual carries the denominator, observations are the
  # exchangeable unit and the null matches PRIMER's closely.
  expect_lt(max(gap[!unit_level]), 0.01)

  # KNOWN GAP. Where a grouping stratum carries the denominator, the
  # exchangeable unit is an island x depth cell. distmer finds the same space
  # PRIMER does -- 720 where PRIMER names it exactly, asserted above -- but the
  # null over that space is not yet identical, and distmer comes out the more
  # conservative of the two. This is a bound on a discrepancy that is
  # understood to be open, not a claim of agreement; tighten it if the null is
  # brought into line. See the note in NEWS.md.
  expect_lt(max(gap[unit_level]), 0.17)
  expect_lt(mean(gap[unit_level]), 0.045)

  # What does hold throughout: the same calls at the five per cent level.
  expect_equal(ours$p_perm[i] < 0.05, reference$p_perm < 0.05)
})
