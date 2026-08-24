# Whole-plot permutation moves intact units, which the block-swap can only do
# between units holding the same number of observations in the same sub-plot
# classes. Real surveys give every unit a distinct count vector, leaving no two
# units exchangeable and every whole-plot term with a single attainable
# arrangement -- p = 1 for any effect size.
#
# `whole_plot` chooses what to do about that.
#   "observation" keeps the intact-unit action at the observation level. Units
#     are weighted by how many observations they carry, which is the estimand
#     the observation-level IS(unit) stratum defines. Exact, but silent when no
#     permutation is possible.
#   "centroid" tests whole-plot terms on equally weighted unit centroids and
#     permutes whole units freely. Unit size no longer restricts exchange, so a
#     p-value exists for any design, but each unit now counts once regardless of
#     effort -- a different estimand, and one that only coincides with
#     "observation" when units are equally replicated.
#   "auto" (default) uses "observation" whenever the intact-unit action can move
#     at least one pair of units, and "centroid" otherwise. The choice is always
#     recorded on the fit and named in each affected test's warnings.
.asym_whole_plot_action <- function(unit_topology, whole_plot) {
  exchangeable <- any(table(unit_topology) > 1L)
  action <- switch(
    whole_plot,
    observation = "observation",
    centroid = "centroid",
    auto = if (exchangeable) "observation" else "centroid"
  )
  list(action = action, exchangeable = exchangeable)
}

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
#' @param denominator How each fixed term's denominator is chosen. `"ems"`
#'   (default) derives the expected mean squares of the numerator and of every
#'   candidate stratum from the actual design, by the trace method, and matches
#'   them with [quasi_denominator()]. When the design is balanced this returns a
#'   single stratum and reproduces `"stratum"` exactly; when it is not, the
#'   correct denominator can be a combination of strata with Satterthwaite
#'   degrees of freedom, and the balanced rule is then wrong. `"stratum"` forces
#'   the balanced rule: whole-plot terms over the unit stratum, sub-plot terms
#'   over the unit-by-sub-plot stratum.
#' @param whole_plot How whole-plot terms are permuted when units are unequally
#'   replicated. `"auto"` (default) keeps the intact-unit action at the
#'   observation level whenever any two units are exchangeable, and otherwise
#'   falls back to equally weighted unit centroids so a p-value exists at all.
#'   `"observation"` and `"centroid"` force either. The two agree exactly for
#'   equally replicated units and target different estimands otherwise; the
#'   choice made is reported on the fit.
#' @param keep_permutations Retain permutation statistics on each test spec.
#' @param max_exact Maximum permutation space allowed for exact enumeration.
#' @param tol Numerical tolerance.
#' @return An object of class `permanova_nested` (shared table/summary shape).
#' @export
permanova_asymmetric_mixed <- function(design, permutations = 999L, seed = NULL,
                                       hypothesis = "I",
                                       denominator = c("ems", "stratum"),
                                       whole_plot = c("auto", "observation", "centroid"),
                                       keep_permutations = FALSE,
                                       max_exact = 100000L, tol = 1e-10) {
  hypothesis <- match.arg(hypothesis, c("I", "II", "III"))
  denominator <- match.arg(denominator)
  whole_plot <- match.arg(whole_plot)
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

  # Trace-derived expected mean squares. The coefficient on each random
  # component is tr(H V) / df with V the component's covariance basis, which is
  # read off the actual design and so stays correct under imbalance. Under
  # balance a numerator's EMS equals one stratum's exactly and the matched
  # denominator collapses to that stratum.
  # The interaction component is orthogonalised against the unit main effect, so
  # that V(unit:subplot) carries only variation the unit term does not already
  # hold. Without it the unit stratum's EMS picks up an interaction coefficient,
  # no random component is attributable to a single stratum, and no exact EMS
  # match exists. This is the standard variance-component parameterisation and
  # reproduces PRIMER's printed EMS, which shows zero there.
  covariance <- list(residual = diag(n))
  covariance[[main_term$label]] <- tcrossprod(design$Z[[main_term$label]])
  interaction_Z <- design$Z[[inter_term$label]]
  interaction_Z <- interaction_Z - Hu %*% interaction_Z
  covariance[[inter_term$label]] <- tcrossprod(interaction_Z)
  ems_of <- function(H, df) {
    vapply(covariance, function(V) trace_product(H, V) / df, numeric(1))
  }
  candidate_ems <- rbind(ems_of(H_is, df_is), ems_of(H_deis, df_deis),
                         ems_of(H_res, df_res))
  rownames(candidate_ems) <- c(paste0("IS(", unit_name, ")"),
                               paste0(paste(subplot_names, collapse = ":"),
                                      " x IS(", unit_name, ")"),
                               "Residual")
  candidate_ms <- c(ms_is, ms_deis, ms_res)
  candidate_df <- c(df_is, df_deis, df_res)

  if (!is.null(seed)) set.seed(seed)
  unit_blocks <- unname(split(seq_len(n), unit))
  subplot_counts <- vapply(unit_blocks, function(indices) {
    as.integer(table(factor(subplot[indices], levels = levels(subplot))))
  }, integer(nlevels(subplot)))
  if (is.null(dim(subplot_counts))) subplot_counts <- matrix(subplot_counts, nrow = nlevels(subplot))
  topology_signature <- apply(subplot_counts, 2L, paste, collapse = ":")
  unit_topology <- factor(topology_signature, levels = unique(topology_signature))
  whole_action <- .asym_whole_plot_action(unit_topology, whole_plot)

  # Unit-level view, used when whole-plot terms are tested on centroids.
  aggregated <- .aggregate_gower_units(design, unit, tol = tol)
  Xu <- X[aggregated$representatives, , drop = FALSE]
  # Unit-level analogue of IS(unit): the residual among unit centroids after
  # every whole-plot term, so all whole-plot tests share one denominator exactly
  # as they do at the observation level.
  unit_full <- projection_matrix(
    Xu[, c(intercept, which(assign %in% which(!involves_subplot))), drop = FALSE],
    tol = tol)
  H_is_unit <- diag(nrow(Xu)) - unit_full$H
  df_is_unit <- nrow(Xu) - unit_full$rank
  ms_is_unit <- trace_hgh(H_is_unit, aggregated$G) / max(df_is_unit, 1L)
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
    quasi <- NULL
    numerator_extra <- 0
    numerator_df_used <- num_df
    if (identical(denominator, "ems")) {
      # Under H0 the fixed term contributes nothing, so its EMS is the residual
      # coefficient plus whatever random components the design leaves in it.
      numerator_ems <- ems_of(num_H, num_df)
      numerator_ems[["residual"]] <- 1
      quasi <- .ems_f_construction(
        numerator_ems, num_ms, num_df,
        candidate_ems, candidate_ms, candidate_df, tol = tol
      )
      if (!is.null(quasi)) {
        denom_ms <- quasi$denominator_ms
        denom_df <- quasi$denominator_df
        num_ms <- quasi$numerator_ms
        numerator_df_used <- quasi$numerator_df
        numerator_extra <- quasi$numerator_residual
        denom_name <- if (quasi$single) {
          quasi$strata_used[[1L]]
        } else {
          paste0("EMS(", paste(quasi$strata_used, collapse = " + "), ")")
        }
        # The permutation action still follows the whole-plot / sub-plot split;
        # only the ratio the statistic is formed from changes.
        if (!quasi$single) denom_H <- NULL
      }
    }
    use_centroid <- identical(kind, "whole") &&
      identical(whole_action$action, "centroid")
    if (use_centroid) {
      # Same hypothesis, evaluated on equally weighted unit centroids so that
      # units of unequal size can be exchanged. The reduced model keeps the
      # other whole-plot terms; sub-plot terms are constant within a unit and
      # drop out of the unit-level design.
      whole_cols <- which(assign %in% which(!involves_subplot))
      adjust_unit <- intersect(adjust_cols, whole_cols)
      reduced_u <- projection_matrix(
        Xu[, c(intercept, adjust_unit), drop = FALSE], tol = tol)
      full_u <- projection_matrix(
        Xu[, c(intercept, adjust_unit, term_cols), drop = FALSE], tol = tol)
      num_df_u <- full_u$rank - reduced_u$rank
      if (num_df_u <= 0L) {
        stop("Whole-plot term `", term, "` is not estimable on unit centroids.",
             call. = FALSE)
      }
      denom_H_u <- H_is_unit
      denom_df_u <- df_is_unit
      if (denom_df_u <= 0L) {
        stop("No unit-centroid residual degrees of freedom remain for `",
             term, "`.", call. = FALSE)
      }
      plan_u <- .permutation_orders(nrow(Xu), permutations, max_exact)
      fl <- .unit_freedman_lane(
        aggregated$G, full_u$H - reduced_u$H, num_df_u, reduced_u$H,
        denom_H_u, denom_df_u, plan_u, tol,
        label = paste0("Permuting term `", term,
                       "` on equally weighted ", unit_name, " centroids"))
      num_ss <- fl$numerator_ss; num_df <- num_df_u; num_ms <- fl$numerator_ms
      num_H <- full_u$H - reduced_u$H
      denom_ms <- fl$denominator_ms; denom_df <- denom_df_u
      denom_name <- paste0(unit_name, "-centroid Residual")
      fl$used <- nrow(plan_u$orders)
      fl$attainable <- plan_u$attainable
      fl$exact <- plan_u$exact && nrow(plan_u$orders) > 0L
    } else {
      fl <- .asym_freedman_lane(G, num_H, num_df, reduced$H, denom_H, denom_df,
                                num_ms / denom_ms, kind, unit, unit_blocks,
                                permutations = permutations,
                                unit_topology = unit_topology,
                                max_exact = max_exact, tol = tol,
                                quasi = if (!is.null(quasi) && !quasi$single) {
                                  list(coefficients = quasi$coefficients,
                                       strata_H = list(H_is, H_deis, H_res),
                                       strata_df = c(df_is, df_deis, df_res))
                                } else NULL,
                                label = paste0(
                                  "Permuting term `", term, "` over ", denom_name,
                                  if (identical(kind, "whole")) {
                                    " (whole units)"
                                  } else " (observations within units)"))
    }
    tested[[term]] <- list(
      term = term,
      distinct_statistics = fl$distinct_statistics,
      minimum_attainable_p = fl$minimum_attainable_p,
      resolution_note = fl$resolution_note,
      numerator_ss = num_ss, numerator_df = num_df,
      numerator_df_used = numerator_df_used, numerator_residual = numerator_extra,
      numerator_ms = num_ms, denominator_ms = denom_ms, denominator_df = denom_df,
      denominator_name = denom_name, statistic = num_ms / denom_ms,
      estimand_note = if (use_centroid) paste0(
        "Units are unequally replicated, so `", term, "` was tested on equally ",
        "weighted ", unit_name, " centroids: intact-unit permutation at the ",
        "observation level had no exchangeable pair. Each ", unit_name,
        " counts once regardless of how many observations it carries, which is ",
        "a different estimand from the observation-level ", is_name,
        " stratum and coincides with it only under equal replication."
      ) else character(),
      whole_plot_action = if (identical(kind, "whole")) {
        if (use_centroid) "centroid" else "observation"
      } else NA_character_,
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
                 distinct_statistics = t$distinct_statistics,
                 minimum_attainable_p = t$minimum_attainable_p,
                 whole_plot_action = t$whole_plot_action,
                 warnings = c(t$estimand_note, t$resolution_note),
                 permutation_statistics = t$permutation_statistics)
    class(spec) <- "permanova_test_spec"; spec
  }), names(tested))

  out <- list(
    whole_plot = whole_action$action,
    whole_plot_exchangeable = whole_action$exchangeable,
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
                                unit_topology, max_exact, tol, label = NULL,
                                quasi = NULL) {
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
  stats_perm <- if (is.null(quasi)) {
    .freedman_lane_statistics(
      G, reduced_H, H_num, num_df, H_denom, denom_df, plan$orders,
      label = label, tol = tol
    )
  } else {
    .freedman_lane_quasi_statistics(
      G, reduced_H, H_num, num_df, quasi$strata_H, quasi$strata_df,
      quasi$coefficients, plan$orders, label = label, tol = tol
    )
  }
  inference <- .permutation_inference(stats_perm, observed, plan$exact, tol)
  p_value <- inference$p_value
  list(p_value = p_value, statistics = stats_perm, used = nrow(plan$orders),
       attainable = plan$attainable, exact = plan$exact && nrow(plan$orders) > 0L,
       distinct_statistics = inference$distinct_statistics,
       minimum_attainable_p = inference$minimum_attainable_p,
       resolution_note = inference$resolution_note)
}
