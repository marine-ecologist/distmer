# =============================================================================
# permanova_contrasts()  -  emmeans-style contrast families for distmer
# -----------------------------------------------------------------------------
# A `type` is a recipe turning factor LEVELS into an L matrix (rows = contrasts,
# cols = levels). Each 1-df contrast is then tested by one of two coherent
# methods:
#
#   method = "isolated" (default)
#       Subset to the involved levels and run the VALIDATED permanova() engine
#       as a two-group test: its own residual, its own restricted permutation.
#       This is what PRIMER's pairwise routine does (manual 1.11). It reuses the
#       omnibus engine directly, so there is NO permutation reconstruction and
#       nothing new to validate.
#
#   method = "pooled"
#       Test the contrast against the FULL-model residual, with numerator
#       SS = trace(H_c G) (deterministically validated: equals the two-group
#       between-SS to machine precision). This answers "is the contrast
#       significant given the whole design's error structure." Its PERMUTATION
#       null must come from the engine's own generator restricted to the
#       contrast's exchangeable units; until that is wired, "pooled" returns the
#       exact statistic and df but REFUSES to fabricate a p-value.
#
# The distinction matters: "pooled" borrows error df (more power, homogeneous-
# dispersion assumption); "isolated" is self-contained (fewer assumptions,
# matches PRIMER). Neither silently uses the full-data permutation space, which
# was the bug in the earlier reconstruction.
#
# INTEGRATION SEAM (section 4): .dm_fixed_factor(), .dm_subset_permanova(),
# .dm_strata_label(), and the not-yet-wired .dm_pooled_contrast_perm().
# =============================================================================

# `%||%`, .stopf, and .warnf are defined once in R/utils.R.

# =============================================================================
# 1 - CONTRAST FAMILIES  -  type -> L matrix (rows = contrasts, cols = levels)
# =============================================================================

#' @keywords internal
#' @noRd
.rows_pairwise <- function(levs) {
  k <- length(levs)
  if (k < 2L) .stopf("pairwise contrasts need >= 2 levels; found %d", k)
  cmb <- utils::combn(k, 2L)
  L   <- matrix(0, ncol(cmb), k, dimnames = list(NULL, levs))
  for (j in seq_len(ncol(cmb))) { L[j, cmb[1, j]] <- 1; L[j, cmb[2, j]] <- -1 }
  rownames(L) <- apply(cmb, 2L, function(ix) paste(levs[ix[1]], "-", levs[ix[2]]))
  L
}

#' @keywords internal
#' @noRd
.rows_trt_vs_ctrl <- function(levs, ref) {
  k <- length(levs); ref <- as.character(ref)
  if (!ref %in% levs) .stopf("ref '%s' is not a level (%s)", ref, paste(levs, collapse = ", "))
  others <- setdiff(levs, ref)
  L <- matrix(0, length(others), k, dimnames = list(NULL, levs))
  for (j in seq_along(others)) { L[j, others[j]] <- 1; L[j, ref] <- -1 }
  rownames(L) <- paste(others, "-", ref)
  L
}

#' @keywords internal
#' @noRd
.rows_poly <- function(levs, scores = NULL) {
  k <- length(levs)
  if (k < 3L) .stopf("polynomial contrasts need >= 3 levels; found %d", k)
  P <- if (is.null(scores)) stats::contr.poly(k) else stats::contr.poly(k, scores = scores)
  L <- t(P); colnames(L) <- levs
  deg <- c("linear", "quadratic", "cubic", "quartic")
  rownames(L) <- ifelse(seq_len(nrow(L)) <= 4L, deg[seq_len(nrow(L))],
                        paste0("degree_", seq_len(nrow(L))))
  L
}

#' @keywords internal
#' @noRd
.rows_contr_sum <- function(levs) {
  k <- length(levs); L <- t(stats::contr.sum(k)); colnames(L) <- levs
  rownames(L) <- paste0(levs[seq_len(k - 1L)], " - mean"); L
}

#' @keywords internal
#' @noRd
.rows_helmert <- function(levs) {
  k <- length(levs); L <- t(stats::contr.helmert(k)); colnames(L) <- levs
  rownames(L) <- paste0("helmert_", seq_len(k - 1L)); L
}

#' @keywords internal
#' @noRd
.rows_custom <- function(levs, custom) {
  k <- length(levs)
  if (is.list(custom) && !is.data.frame(custom)) {
    L <- do.call(rbind, custom)
    if (is.null(rownames(L))) rownames(L) <- names(custom)
  } else L <- as.matrix(custom)
  if (ncol(L) != k) .stopf("custom contrast has %d cols but factor has %d levels", ncol(L), k)
  colnames(L) <- levs
  if (any(abs(rowSums(L)) > 1e-8)) .warnf("some custom rows do not sum to zero")
  if (is.null(rownames(L))) rownames(L) <- paste0("contrast_", seq_len(nrow(L)))
  L
}

#' Build a contrast (L) matrix from a named family
#' @param type One of "pairwise","trt.vs.ctrl","poly","contr.sum","helmert",
#'   "custom", or a matrix / named list of contrast rows.
#' @param levels Character vector of factor levels, in model order.
#' @param ref Reference level for "trt.vs.ctrl".
#' @param scores Optional numeric spacing for "poly".
#' @param custom Matrix or named list of rows when type = "custom".
#' @return Numeric matrix: rows = contrasts (named), columns = levels (named).
#' @export
contrast_matrix <- function(type, levels, ref = NULL, scores = NULL, custom = NULL) {
  levs <- as.character(levels)
  if (is.matrix(type) || (is.list(type) && !is.data.frame(type))) return(.rows_custom(levs, type))
  switch(
    match.arg(type, c("pairwise", "trt.vs.ctrl", "poly", "contr.sum", "helmert", "custom")),
    pairwise    = .rows_pairwise(levs),
    trt.vs.ctrl = .rows_trt_vs_ctrl(levs, ref %||% levs[1L]),
    poly        = .rows_poly(levs, scores = scores),
    contr.sum   = .rows_contr_sum(levs),
    helmert     = .rows_helmert(levs),
    custom      = .rows_custom(levs, custom %||% .stopf("type='custom' needs `custom`"))
  )
}

# which levels does a contrast row actually involve (nonzero weight)? ----------
#' @keywords internal
#' @noRd
.contrast_levels <- function(L_row, levs) levs[abs(L_row) > 1e-8]


# =============================================================================
# 2 - MAIN ENTRY  -  permanova_contrasts()
# =============================================================================

#' Contrasts for a fitted PERMANOVA
#'
#' Planned or pairwise contrasts on a fixed factor of a fitted [permanova()]
#' model. Contrast families mirror \pkg{emmeans}. Two testing methods are
#' offered (see `method`), differing in denominator and permutation scope.
#'
#' @param fit A `permanova` object.
#' @param type Contrast family: "pairwise", "trt.vs.ctrl", "poly",
#'   "contr.sum", "helmert", or a contrast matrix / named list.
#' @param method Testing method. `"isolated"` (default) runs each contrast as a
#'   two-group [permanova()] on the involved levels only - its own residual and
#'   restricted permutation; this matches PRIMER's pairwise routine and reuses
#'   the validated engine. `"pooled"` tests against the full-model residual with
#'   numerator SS = trace(H_c G); it borrows error df (more power) but assumes
#'   homogeneous dispersion, and its permutation null requires the engine's own
#'   generator (see Details).
#' @param factor Name of the fixed factor. Defaults to the tested factor.
#' @param ref Reference level for "trt.vs.ctrl".
#' @param scores Numeric spacing for "poly".
#' @param custom Matrix / named list for type = "custom".
#' @param adjust Multiplicity method for [stats::p.adjust()] (default "holm").
#' @param permutations Number of permutations; defaults to `fit`'s value.
#' @param seed RNG seed; defaults to `fit`'s value.
#'
#' @details
#' `"isolated"` is exact and fully validated: it is literally [permanova()] on
#' the subset, so its statistic, denominator and permutation match a two-group
#' analysis of those levels. `"pooled"` uses a numerator SS that is validated
#' deterministically (it equals the between-group SS to machine precision), but
#' its permutation null must be produced by the engine's Freedman-Lane generator
#' restricted to the contrast's exchangeable units. Until that generator is
#' wired in, `"pooled"` returns the statistic, df and denominator but leaves the
#' p-value `NA` rather than computing it over the wrong (full-data) permutation
#' space.
#'
#' @return An object of class `permanova_contrasts`.
#' @seealso [permanova()], [contrast_matrix()]
#' @export
permanova_contrasts <- function(fit,
                                type         = "pairwise",
                                method       = c("isolated", "pooled"),
                                factor       = NULL,
                                ref          = NULL,
                                scores       = NULL,
                                custom       = NULL,
                                adjust       = "holm",
                                permutations = NULL,
                                seed         = NULL) {

  if (!inherits(fit, "permanova")) .stopf("`fit` must be a permanova() object")
  method       <- match.arg(method)
  adjust       <- match.arg(adjust, stats::p.adjust.methods)
  permutations <- permutations %||% fit$permutations %||% fit$mc_permutations %||% 999L
  seed         <- seed         %||% fit$seed

  f <- .dm_fixed_factor(fit, factor)
  if (length(f$levels) < 2L) .stopf("factor '%s' has < 2 levels", f$name)

  if (identical(type, "poly") && is.null(scores)) {
    num <- suppressWarnings(as.numeric(f$levels))
    if (!any(is.na(num)) && length(unique(num)) == length(num)) {
      scores <- num
      if (!f$is_ordered)
        .warnf("factor '%s' not ordered; using level values as poly scores", f$name)
    }
  }

  L <- contrast_matrix(if (!is.null(custom)) custom else type,
                       levels = f$levels, ref = ref, scores = scores, custom = custom)

  rows <- lapply(seq_len(nrow(L)), function(i) {
    r <- if (method == "isolated") {
      .contrast_isolated(fit, f, L[i, ], permutations, seed)
    } else {
      .contrast_pooled(fit, f, L[i, ], permutations, seed)
    }
    data.frame(
      contrast     = rownames(L)[i],
      pseudo_F     = r$pseudo_F,
      df_num       = r$df_num       %||% 1L,
      df_den       = r$df_den       %||% NA_real_,
      denominator  = r$denom_term   %||% NA_character_,
      p_perm       = r$p_perm,
      unique_perms = r$unique_perms %||% NA_real_,
      estimable    = r$estimable    %||% TRUE,
      stringsAsFactors = FALSE
    )
  })
  tab <- do.call(rbind, rows)

  est <- tab$estimable & !is.na(tab$p_perm)
  tab$p_adj <- NA_real_
  if (any(est)) tab$p_adj[est] <- stats::p.adjust(tab$p_perm[est], method = adjust)

  .coarse_warn(tab, permutations)
  if (method == "pooled" && all(is.na(tab$p_perm)))
    .warnf("method='pooled': permutation p-values not computed (engine generator not wired). Statistic and df are valid; use method='isolated' for tested p-values.")

  structure(list(
    table = tab, L = L, factor = f$name, levels = f$levels,
    type = if (is.character(type)) type else "custom",
    method = method, adjust = adjust, permutations = permutations,
    strata = .dm_strata_label(fit), seed = seed
  ), class = "permanova_contrasts")
}

#' @keywords internal
#' @noRd
.coarse_warn <- function(tab, permutations) {
  up <- suppressWarnings(as.numeric(tab$unique_perms))
  if (all(is.na(up))) return(invisible())
  requested <- if (is.numeric(permutations)) permutations + 1L else Inf
  floorp <- 1 / pmin(up, requested)
  hit    <- which(!is.na(floorp) & floorp > 0.05)
  if (length(hit))
    .warnf("contrast(s) %s have a coarse permutation space; smallest attainable p ~ %.3g",
           paste(tab$contrast[hit], collapse = "; "), max(floorp[hit], na.rm = TRUE))
  invisible()
}


# =============================================================================
# 3 - THE TWO METHODS
# =============================================================================

# isolated: subset to the involved levels, run the validated engine as a
# two-group test. Exact, PRIMER-matching, no reconstruction.
#' @keywords internal
#' @noRd
.contrast_isolated <- function(fit, f, L_row, permutations, seed) {
  involved <- .contrast_levels(L_row, f$levels)
  if (length(involved) != 2L) {
    # multi-level contrasts (poly, contr.sum) have no clean two-group subset;
    # isolated is defined only for 2-level contrasts. Fall back to pooled-stat
    # with NA p and a note carried via denom_term.
    st <- .contrast_pooled(fit, f, L_row, permutations, seed)
    st$denom_term <- "Residual (pooled; isolated undefined for >2-level contrast)"
    return(st)
  }
  sub_fit <- .dm_subset_permanova(fit, f$name, involved, permutations, seed)
  sub <- sub_fit$table
  # sub is the two-group permanova table; pull its factor row
  frow <- sub[sub$Source == f$name | grepl(paste0("^", f$name), sub$Source), , drop = FALSE][1, ]
  sub_tests <- .permanova_tests(sub_fit$fit)
  spec <- sub_tests[[f$name]]
  if (is.null(spec) && length(sub_tests)) spec <- sub_tests[[1L]]
  pcol <- if ("P(perm)" %in% names(sub)) "P(perm)" else "p_perm"
  ucol <- if ("Unique perms" %in% names(sub)) "Unique perms" else "unique_perms"
  list(
    pseudo_F     = frow[["Pseudo-F"]] %||% frow[["pseudo_F"]],
    df_num       = frow[["df"]]       %||% frow[["df_num"]] %||% 1L,
    df_den       = spec$denominator$df %||% NA_real_,
    denom_term   = spec$denominator$expression %||% "two-group denominator",
    p_perm       = frow[[pcol]],
    unique_perms = frow[[ucol]],
    estimable    = TRUE
  )
}

# pooled: numerator SS = trace(H_c G) tested against the full-model residual.
# Statistic + df are exact; p-value requires the engine's restricted generator.
#' @keywords internal
#' @noRd
.contrast_pooled <- function(fit, f, L_row, permutations, seed) {
  design <- fit$design
  G      <- design$G
  tinfo  <- design$projections$terms[[f$name]] %||%
    .stopf("no per-term projection for '%s'", f$name)

  # rank-1 contrast regressor within the term space, residualised by H_reduced
  fcol <- design$data[[f$name]]; if (!is.factor(fcol)) fcol <- factor(fcol)
  Ind  <- stats::model.matrix(~ 0 + fcol); colnames(Ind) <- levels(fcol)
  w    <- L_row[levels(fcol)]
  xc_raw <- Ind %*% w
  Hr   <- tinfo$H_reduced
  xc   <- xc_raw - Hr %*% xc_raw

  scale0 <- sqrt(sum(xc_raw^2)); normc <- sqrt(sum(xc^2))
  estimable <- is.finite(normc) && scale0 > 0 && (normc / scale0) > 1e-8
  if (!estimable)
    return(list(pseudo_F = NA_real_, p_perm = NA_real_, df_num = 1L, df_den = NA_real_,
                denom_term = "Residual", unique_perms = NA_real_, estimable = FALSE))

  denomSS <- as.numeric(crossprod(xc))
  SS_c    <- if (denomSS > 0) as.numeric(crossprod(xc, G %*% xc)) / denomSS else 0
  tab     <- fit$table
  resid   <- tab[grepl("Residual", tab$Source), , drop = FALSE][1, ]
  MS_res  <- resid[["MS"]]; df_res <- resid[["df"]]
  F_obs   <- (SS_c / 1) / MS_res

  # permutation null: MUST come from the engine's restricted generator.
  perm <- .dm_pooled_contrast_perm(fit, xc, permutations = permutations, seed = seed)

  list(
    pseudo_F     = F_obs,
    df_num       = 1L,
    df_den       = df_res,
    denom_term   = "Residual (full model, pooled)",
    p_perm       = perm$p_perm,          # NA until the generator is wired
    unique_perms = perm$unique_perms,    # NA until the generator is wired
    estimable    = TRUE
  )
}


# =============================================================================
# 4 - PRINT / SUMMARY / as.data.frame
# =============================================================================

#' @param x A `permanova_contrasts` object.
#' @param digits Significant digits.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @noRd
#' @exportS3Method base::print
print.permanova_contrasts <- function(x, digits = 3, ...) {
  cat(sprintf("Permutational contrasts on factor '%s'\n", x$factor))
  cat(sprintf("  family      : %s  (%d contrasts)\n", x$type, nrow(x$table)))
  cat(sprintf("  method      : %s\n", x$method))
  cat(sprintf("  levels      : %s\n", paste(x$levels, collapse = ", ")))
  cat(sprintf("  exchangeable: %s\n", x$strata %||% "unrestricted"))
  cat(sprintf("  permutations: %s   p-adjust: %s\n\n", x$permutations, x$adjust))
  tb <- x$table
  show <- data.frame(
    contrast    = tb$contrast,
    pseudo_F    = round(tb$pseudo_F, digits),
    df          = ifelse(is.na(tb$df_den), as.character(tb$df_num),
                         paste(tb$df_num, tb$df_den, sep = ", ")),
    denominator = tb$denominator,
    p_perm      = ifelse(is.na(tb$p_perm), "-",
                         format.pval(tb$p_perm, digits = digits, eps = 1e-4)),
    p_adj       = ifelse(is.na(tb$p_adj), "-",
                         format.pval(tb$p_adj, digits = digits, eps = 1e-4)),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  print(show, row.names = FALSE, right = FALSE)
  if (x$method == "pooled" && all(is.na(tb$p_perm)))
    cat("\n  note: pooled p-values require the engine permutation generator;",
        "statistic and df are exact. Use method='isolated' for tested p-values.\n")
  invisible(x)
}

#' @param object A `permanova_contrasts` object.
#' @param ... Unused.
#' @return The results `data.frame`.
#' @noRd
#' @exportS3Method base::summary
summary.permanova_contrasts <- function(object, ...) object$table

#' @param x A `permanova_contrasts` object.
#' @param row.names,optional,... Standard `as.data.frame` args (unused).
#' @return The results `data.frame`.
#' @noRd
#' @exportS3Method base::as.data.frame
as.data.frame.permanova_contrasts <- function(x, row.names = NULL, optional = FALSE, ...) x$table


# =============================================================================
# 5 - INTEGRATION SEAM  -  bind to distmer internals
# =============================================================================

# resolve fixed factor + levels from the compiled design
#' @keywords internal
#' @noRd
.dm_fixed_factor <- function(fit, factor = NULL) {
  design <- fit$design %||% .stopf("fit has no $design")
  dat    <- design$data %||% .stopf("fit$design has no $data")
  fx     <- design$fixed_effects
  nm <- factor %||% {
    if (length(fx) < 1L) .stopf("no fixed effects; pass `factor =`")
    if (length(fx) > 1L)
      .warnf("multiple fixed effects (%s); contrasting '%s'. Use `factor=`.",
             paste(fx, collapse = ", "), fx[[1L]])
    fx[[1L]]
  }
  if (!nm %in% names(dat)) .stopf("factor '%s' not in fit$design$data", nm)
  col <- dat[[nm]]; if (!is.factor(col)) col <- factor(col)
  list(name = nm, levels = levels(col), values = col, is_ordered = is.ordered(col))
}

# isolated engine: subset to the involved levels and RE-RUN permanova().
# This is the validated path - it calls the real engine, not a reconstruction.
#' @keywords internal
#' @noRd
.dm_subset_permanova <- function(fit, factor_name, levels_keep, permutations, seed) {
  # REAL implementation (uncomment/adapt to your call signature):
  #
  # design <- fit$design
  # keep   <- design$data[[factor_name]] %in% levels_keep
  # d_full <- as.dist(design$distance)                 # or stored resemblance
  # d_sub  <- as.dist(as.matrix(d_full)[keep, keep])
  # env    <- droplevels(design$data[keep, , drop = FALSE])
  # form   <- stats::reformulate(factor_name, response = "d_sub")
  # permanova(form, data = env, permutations = permutations, seed = seed)$table
  #
  .stopf(paste0(".dm_subset_permanova(): wire to a subset re-run of permanova() ",
                "(subset the resemblance + data to `levels_keep`, refit). ",
                "This reuses the validated engine directly."))
}

# pooled permutation: MUST be the engine's Freedman-Lane generator, restricted
# to the contrast's exchangeable units and fed the rank-1 regressor `xc`.
# Returns NA until wired, so pooled never fabricates a p over the wrong space.
#' @keywords internal
#' @noRd
.dm_pooled_contrast_perm <- function(fit, xc, permutations, seed = NULL) {
  # REAL implementation: call whatever produced fit$fit$permutation, passing the
  # rank-1 numerator projection built from `xc`, so the null uses the SAME
  # restricted exchangeable units as the omnibus term. Return:
  #   list(p_perm = <value>, unique_perms = <count for THIS hypothesis>)
  list(p_perm = NA_real_, unique_perms = NA_real_)   # honest NA, not a wrong p
}

#' @keywords internal
#' @noRd
.dm_strata_label <- function(fit) {
  fit$fit$permutation$algorithm %||% fit$fit$supported_design %||%
    "full-model residual"
}


# =============================================================================
# .dm_subset_permanova()  -  wired to the validated permanova() engine
# -----------------------------------------------------------------------------
# Replaces the stub in permanova_contrasts.R section 5. For a two-level
# contrast, subset the resemblance + design data to the involved levels and
# re-run permanova() as a two-group test. This reuses the engine you have
# already validated against PRIMER (A1/A2), so the contrast inherits the exact
# statistic, denominator, and permutation scheme - no reconstruction.
# =============================================================================

#' @keywords internal
#' @noRd
.dm_subset_permanova <- function(fit, factor_name, levels_keep,
                                 permutations, seed) {

  design <- fit$design %||% .stopf("fit has no $design")

  # 1. the full resemblance the fit was built on
  #    (fit$design$distance is the n x n matrix; fall back to $G-derived if absent)
  D_full <- design$distance %||% .stopf("fit$design has no $distance")
  Dm     <- as.matrix(D_full)

  # 2. rows/cols to keep: samples whose factor level is in levels_keep
  fac  <- design$data[[factor_name]]
  if (!is.factor(fac)) fac <- factor(fac)
  keep <- as.character(fac) %in% as.character(levels_keep)
  if (sum(keep) < 2L)
    .stopf("subset for %s has < 2 samples", paste(levels_keep, collapse = " vs "))

  # 3. subset the resemblance and the design data, drop unused factor levels
  d_sub   <- stats::as.dist(Dm[keep, keep, drop = FALSE])
  env_sub <- droplevels(design$data[keep, , drop = FALSE])

  args <- list(formula = design$formula, data = env_sub, distance = d_sub,
               permutations = permutations, seed = seed)
  if (!is.null(fit$type) && length(fit$type) == 1L && !is.na(fit$type)) {
    args$type <- fit$type
  }
  do.call(permanova, args)
}

# =============================================================================
# 6 - USAGE
# =============================================================================
# fit <- permanova(pyr.bray ~ dose, data = pyr.env, permutations = 999)
#
# # PRIMER-matching pairwise (default): each pair as its own two-group test
# permanova_contrasts(fit, type = "pairwise")                  # method="isolated"
# permanova_contrasts(fit, type = "trt.vs.ctrl", ref = "0")
#
# # pooled against the full-model residual (more power); p requires the engine
# # generator - statistic + df are exact, p is NA until wired
# permanova_contrasts(fit, type = "pairwise", method = "pooled")
