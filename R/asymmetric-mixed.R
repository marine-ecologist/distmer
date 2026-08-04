#' Asymmetric nested mixed PERMANOVA (unit + unit-by-subplot random strata)
#'
#' Reproduces PRIMER's asymmetric nested mixed design with two random terms: a
#' unit nested in the crossed whole-plot factors, `(1 | unit)`, and the unit by
#' sub-plot-factor interaction, `(1 | unit:subplot)`. This gives three error
#' strata, all at the observation level (as in PRIMER, which permutes residuals
#' under a reduced model):
#'
#' * **IS(unit)** = between-unit variation after the whole-plot fixed terms --
#'   the denominator for the whole-plot fixed effects (those not involving the
#'   sub-plot factor);
#' * **subplot x IS(unit)** = the unit-by-sub-plot interaction -- the denominator
#'   for the sub-plot fixed effects (those involving the sub-plot factor);
#' * **Residual** = within unit-by-sub-plot cells -- the denominator for both
#'   random strata.
#'
#' Degrees of freedom are derived from the actual (possibly incomplete) design,
#' so units missing a sub-plot level are handled. For the Caribbean benthos
#' example this is `bray ~ DEPTH * ECOREGION * POP_STATUS + (1 | ISLAND) +
#' (1 | ISLAND:DEPTH)`, giving `IS(EC x PS)` and `DE x IS(EC x PS)`.
#'
#' @param design A compiled Design Object with two random terms: one `(1 | unit)`
#'   main term and one `(1 | unit:subplot)` interaction term.
#' @param permutations A non-negative integer or `"all"`.
#' @param seed Optional random seed.
#' @param hypothesis Sum-of-squares partition for the fixed terms (`"I"`, `"II"`,
#'   `"III"`).
#' @param keep_permutations Retain permutation statistics on each test spec.
#' @param max_exact Maximum permutation space allowed for exact enumeration.
#' @param tol Numerical tolerance.
#' @return An object of class `permanova_nested` (shared table/summary shape).
#' @export
permanova_asymmetric_mixed <- function(design, permutations = 999L, seed = NULL,
                                       hypothesis = "I", keep_permutations = FALSE,
                                       max_exact = 100000L, tol = 1e-10) {
  hypothesis <- match.arg(hypothesis, c("I", "II", "III"))
  if (!inherits(design, "permanova_design") || is.null(design$G)) {
    stop("`design` must be a compiled Design Object containing `distance`.",
         call. = FALSE)
  }
  origins <- vapply(design$random_terms, `[[`, character(1), "origin")
  if (length(origins) != 2L || sum(origins == "main") != 1L ||
      sum(origins == "interaction") != 1L) {
    stop("Asymmetric mixed inference requires one `(1 | unit)` term and one ",
         "`(1 | unit:subplot)` interaction term.", call. = FALSE)
  }
  main_term <- design$random_terms[[which(origins == "main")]]
  inter_term <- design$random_terms[[which(origins == "interaction")]]
  unit_name <- main_term$variables
  subplot_names <- setdiff(inter_term$variables, unit_name)
  if (length(unit_name) != 1L || !length(subplot_names) ||
      !all(unit_name %in% inter_term$variables)) {
    stop("The interaction random term must be the unit crossed with one or more ",
         "sub-plot factors, `(1 | unit:subplot)`.", call. = FALSE)
  }
  if (identical(hypothesis, "III")) .validate_type3_contrasts(design, tol = tol)

  n <- design$n
  G <- design$G
  unit <- droplevels(interaction(design$data[unit_name], drop = TRUE))
  cell <- droplevels(interaction(design$data[c(unit_name, subplot_names)], drop = TRUE))
  subplot <- droplevels(interaction(design$data[subplot_names], drop = TRUE))

  Hu <- projection_matrix(stats::model.matrix(~ 0 + unit), tol = tol)$H
  Hc <- projection_matrix(stats::model.matrix(~ 0 + cell), tol = tol)$H
  rank_unit <- as.integer(round(sum(diag(Hu))))
  rank_cell <- as.integer(round(sum(diag(Hc))))

  X <- design$X
  assign <- attr(X, "assign")
  labels <- attr(stats::terms(design$fixed_formula), "term.labels")
  intercept <- which(assign == 0L)
  term_variables <- lapply(labels, function(term) strsplit(term, ":", fixed = TRUE)[[1L]])
  # A fixed term is a sub-plot term when it involves any sub-plot factor,
  # otherwise a whole-plot term.
  involves_subplot <- vapply(term_variables, function(v) any(v %in% subplot_names),
                             logical(1))
  between_terms <- labels[!involves_subplot]
  within_terms <- labels[involves_subplot]
  between_columns <- which(assign %in% which(!involves_subplot))
  within_columns <- which(assign %in% which(involves_subplot))

  # --- the three random error strata (observation level) --------------------
  is_reduced <- projection_matrix(X[, c(intercept, between_columns), drop = FALSE], tol = tol)
  H_is <- Hu - is_reduced$H                                  # IS(unit)
  df_is <- rank_unit - is_reduced$rank
  deis_reduced <- projection_matrix(
    cbind(stats::model.matrix(~ 0 + unit), X[, within_columns, drop = FALSE]), tol = tol
  )
  H_deis <- Hc - deis_reduced$H                              # subplot x IS(unit)
  df_deis <- rank_cell - deis_reduced$rank
  H_res <- diag(n) - Hc                                      # residual
  df_res <- n - rank_cell
  for (nm in c("df_is", "df_deis", "df_res")) {
    if (get(nm) <= 0L) stop("Stratum `", nm, "` has no degrees of freedom; the ",
                            "design does not identify the three-stratum model.",
                            call. = FALSE)
  }
  ss_is <- trace_hgh(H_is, G); ms_is <- ss_is / df_is
  ss_deis <- trace_hgh(H_deis, G); ms_deis <- ss_deis / df_deis
  ss_res <- trace_hgh(H_res, G); ms_res <- ss_res / df_res

  if (!is.null(seed)) set.seed(seed)
  unit_blocks <- unname(split(seq_len(n), unit))
  subplot_counts <- vapply(unit_blocks, function(indices) {
    as.integer(table(factor(subplot[indices], levels = levels(subplot))))
  }, integer(nlevels(subplot)))
  if (is.null(dim(subplot_counts))) subplot_counts <- matrix(subplot_counts, nrow = nlevels(subplot))
  topology_signature <- apply(subplot_counts, 2L, paste, collapse = ":")
  unit_topology <- factor(topology_signature, levels = unique(topology_signature))
  tested <- list()

  # --- fixed-term tests: the numerator SS is standard full-model sequential
  # (Type I/II/III over every preceding fixed term); the whole-plot vs sub-plot
  # classification only selects the denominator stratum.
  is_name <- paste0("IS(", unit_name, ")")
  deis_name <- paste0(paste(subplot_names, collapse = ":"), " x IS(", unit_name, ")")
  for (index in seq_along(labels)) {
    term <- labels[[index]]
    adjust <- labels[.ss_type_adjust(index, labels, hypothesis)]
    adjust_cols <- which(assign %in% match(adjust, labels))
    term_cols <- which(assign == index)
    reduced <- projection_matrix(X[, c(intercept, adjust_cols), drop = FALSE], tol = tol)
    full <- projection_matrix(X[, c(intercept, adjust_cols, term_cols), drop = FALSE], tol = tol)
    num_H <- full$H - reduced$H
    num_df <- full$rank - reduced$rank
    if (num_df <= 0L) stop("Term `", term, "` is not estimable for Type ",
                           hypothesis, " inference.", call. = FALSE)
    num_ss <- trace_hgh(num_H, G); num_ms <- num_ss / num_df
    if (involves_subplot[[index]]) {
      denom_H <- H_deis; denom_df <- df_deis; denom_ms <- ms_deis
      denom_name <- deis_name; kind <- "within"
    } else {
      denom_H <- H_is; denom_df <- df_is; denom_ms <- ms_is
      denom_name <- is_name; kind <- "whole"
    }
    fl <- .asym_freedman_lane(G, num_H, num_df, reduced$H, denom_H, denom_df,
                              num_ms / denom_ms, kind, unit, unit_blocks,
                              permutations = permutations,
                              unit_topology = unit_topology,
                              max_exact = max_exact, tol = tol)
    tested[[term]] <- list(
      term = term, numerator_ss = num_ss, numerator_df = num_df,
      numerator_ms = num_ms, denominator_ms = denom_ms, denominator_df = denom_df,
      denominator_name = denom_name, statistic = num_ms / denom_ms,
      p_value = fl$p_value, permutations_used = fl$used,
      attainable_permutations = fl$attainable,
      exact_permutation = fl$exact,
      permutation_statistics = if (keep_permutations) fl$statistics else NULL)
  }

  # --- random strata tested over the residual (no permutation p by default) --
  strata <- list(
    list(term = paste0("IS(", unit_name, ")"), ss = ss_is, df = df_is, ms = ms_is,
         f = ms_is / ms_res),
    list(term = paste0(paste(subplot_names, collapse = ":"), " x IS(", unit_name, ")"),
         ss = ss_deis, df = df_deis, ms = ms_deis, f = ms_deis / ms_res)
  )

  fixed_rows <- do.call(rbind, lapply(tested, function(t) data.frame(
    term = t$term, df = t$numerator_df, ss = t$numerator_ss, ms = t$numerator_ms,
    denominator = paste0("MS[", t$denominator_name, "]"),
    denominator_df = t$denominator_df, pseudo_f = t$statistic, p_perm = t$p_value,
    stringsAsFactors = FALSE)))
  stratum_rows <- do.call(rbind, lapply(strata, function(s) data.frame(
    term = s$term, df = s$df, ss = s$ss, ms = s$ms, denominator = "MS[Residual]",
    denominator_df = df_res, pseudo_f = s$f, p_perm = NA_real_,
    stringsAsFactors = FALSE)))
  table <- rbind(fixed_rows, stratum_rows, data.frame(
    term = "Residual", df = df_res, ss = ss_res, ms = ms_res,
    denominator = NA_character_, denominator_df = NA_real_, pseudo_f = NA_real_,
    p_perm = NA_real_))
  rownames(table) <- NULL

  specs <- stats::setNames(lapply(tested, function(t) {
    spec <- list(term = t$term,
                 numerator = list(ss = t$numerator_ss, df = t$numerator_df, ms = t$numerator_ms),
                 denominator = list(expression = paste0("MS[", t$denominator_name, "]"),
                                    ms = t$denominator_ms, df = t$denominator_df),
                 statistic = t$statistic, p_value = t$p_value,
                 permutations_used = t$permutations_used,
                 attainable_permutations = t$attainable_permutations,
                 exact_permutation = t$exact_permutation,
                 permutation_statistics = t$permutation_statistics)
    class(spec) <- "permanova_test_spec"; spec
  }), names(tested))

  out <- list(
    call = match.call(), design = design, table = table, tests = specs,
    hierarchy = c(unit_name, subplot_names),
    supported_design = paste(
      "asymmetric nested mixed design: whole-plot fixed effects over the unit",
      "stratum, sub-plot fixed effects over the unit-by-sub-plot stratum, both",
      "random strata over the residual"),
    permutation = list(
      algorithm = "Freedman-Lane residual permutation at the whole-unit and within-unit level",
      used = max(vapply(specs, `[[`, numeric(1), "permutations_used"), 0),
      exact = all(vapply(specs, `[[`, logical(1), "exact_permutation")), seed = seed)
  )
  class(out) <- "permanova_nested"
  out
}

# Freedman-Lane for the asymmetric mixed engine. Whole-plot terms permute whole
# units; sub-plot terms permute observations within units. Both recompute the
# numerator and its stratum denominator under the permutation.
.asym_freedman_lane <- function(G, H_num, num_df, reduced_H, H_denom, denom_df,
                                observed, kind, unit, unit_blocks, permutations,
                                unit_topology, max_exact, tol, label = NULL) {
  n <- nrow(G)
  plan <- if (identical(kind, "whole")) {
    unit_plan <- .stratified_permutation_plan(unit_topology, permutations, max_exact)
    orders <- if (nrow(unit_plan$orders)) {
      t(vapply(seq_len(nrow(unit_plan$orders)), function(i) {
        .block_permutation_index(unit_blocks, unit_plan$orders[i, ])
      }, integer(n)))
    } else matrix(integer(), nrow = 0L, ncol = n)
    list(orders = orders, exact = unit_plan$exact, attainable = unit_plan$attainable)
  } else {
    within_plan <- .stratified_permutation_plan(unit, permutations, max_exact)
    list(orders = within_plan$orders, exact = within_plan$exact,
         attainable = within_plan$attainable)
  }
  stats_perm <- .freedman_lane_statistics(
    G, reduced_H, H_num, num_df, H_denom, denom_df, plan$orders,
    label = label, tol = tol
  )
  p_value <- if (length(stats_perm)) {
    extreme <- sum(stats_perm >= observed - tol, na.rm = TRUE)
    if (plan$exact) extreme / length(stats_perm) else (extreme + 1) / (length(stats_perm) + 1)
  } else NA_real_
  list(p_value = p_value, statistics = stats_perm, used = nrow(plan$orders),
       attainable = plan$attainable, exact = plan$exact && nrow(plan$orders) > 0L)
}
