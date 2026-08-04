# Classify each fixed term relative to a single random unit, for the split-plot
# engines. A pure-factor term is whole_plot (constant within every unit),
# within_unit (fully crossed within every unit, or an interaction with any
# within-unit component), or partial (neither -> confounded with the unit). A
# covariate term is within_unit when it varies within units, else a between-unit
# (whole-plot) covariate. Returns the classification and the column index sets.
.classify_split_plot_terms <- function(design, group, covariates, H_within,
                                       tol = 1e-10) {
  X <- design$X
  assign <- attr(X, "assign")
  labels <- attr(stats::terms(design$fixed_formula), "term.labels")
  term_variables <- lapply(labels, function(term) strsplit(term, ":", fixed = TRUE)[[1L]])
  term_has_covariate <- vapply(term_variables, function(v) any(v %in% covariates),
                               logical(1))
  # Classify each term by how much of its rank lives within the units, so a
  # design need not be fully crossed within every unit -- only estimable:
  #   * within-unit rank 0            -> whole_plot (tested over the unit stratum)
  #   * within-unit rank == term rank -> within_unit (tested over the residual),
  #     e.g. a sub-plot factor that is unbalanced or missing from some units but
  #     whose whole effect is still recovered within units;
  #   * 0 < within-unit rank < term rank -> partial: some contrast is confounded
  #     with the unit and no single stratum identifies it, so it is refused.
  term_class <- vapply(seq_along(labels), function(j) {
    cols <- which(assign == j)
    if (!length(cols)) return("whole_plot")
    Xc <- X[, cols, drop = FALSE]
    term_rank <- projection_matrix(Xc, tol = tol)$rank
    within_rank <- .projection_in_subspace(H_within, Xc, tol = tol)$rank
    if (within_rank == 0L) "whole_plot"
    else if (within_rank >= term_rank) "within_unit"
    else "partial"
  }, character(1))
  is_within <- term_class == "within_unit"
  list(
    labels = labels, assign = assign, term_variables = term_variables,
    term_has_covariate = term_has_covariate, term_class = term_class,
    is_within = is_within, partial_terms = labels[term_class == "partial"],
    between_terms = labels[!is_within], within_terms = labels[is_within],
    between_columns = which(assign %in% which(!is_within)),
    within_columns = which(assign %in% which(is_within))
  )
}

#' Unbalanced multi-factor split-plot PERMANOVA
#'
#' Tests several categorical factors in an unequally replicated `(1 | unit)`
#' random-intercept design. Whole-plot factors (constant within each unit) are
#' tested on **equally weighted unit centroids** over the between-unit stratum --
#' the unit-weighted estimand of [permanova_unequal_random_intercept()],
#' generalised to any number of fixed factors and sums-of-squares type -- with
#' size-restricted whole-unit permutation. Within-unit (sub-plot) factors, and
#' any between-by-within interaction, are tested over the within-unit residual
#' with within-unit residual permutation. This is the distance-based analogue of
#' PRIMER's asymmetric nested design (for example `IS(EC x PS)` as the whole-plot
#' stratum and `DE x IS(EC x PS)` as the sub-plot residual) when islands are
#' unequally sampled.
#'
#' @param design A compiled Design Object with one `(1 | unit)` random term and
#'   one or more categorical fixed factors.
#' @param permutations A non-negative integer or `"all"`.
#' @param seed Optional random seed.
#' @param hypothesis Sum-of-squares partition for the fixed terms (`"I"`, `"II"`,
#'   `"III"`), applied within the whole-plot and within-unit blocks separately.
#' @param keep_permutations Retain permutation statistics on each test spec.
#' @param max_exact Maximum permutation space allowed for exact enumeration.
#' @param tol Numerical tolerance.
#' @return An object of class `permanova_nested` (shared table/summary shape).
#' @export
permanova_unbalanced_splitplot <- function(design, permutations = 999L,
                                           seed = NULL, hypothesis = "I",
                                           keep_permutations = FALSE,
                                           max_exact = 100000L, tol = 1e-10) {
  hypothesis <- match.arg(hypothesis, c("I", "II", "III"))
  if (!inherits(design, "permanova_design") || is.null(design$G)) {
    stop("`design` must be a compiled Design Object containing `distance`.",
         call. = FALSE)
  }
  if (length(design$random_terms) != 1L ||
      design$random_terms[[1L]]$origin != "main") {
    stop("Unbalanced split-plot inference requires one `(1 | unit)` term.",
         call. = FALSE)
  }
  if (identical(hypothesis, "III")) .validate_type3_contrasts(design, tol = tol)
  split <- .nested_covariate_split(design)

  group_label <- design$random_terms[[1L]]$label
  group <- droplevels(factor(design$groups[[group_label]]))
  n <- design$n
  G <- design$G
  n_units <- nlevels(group)

  Z <- stats::model.matrix(~ 0 + group)
  H_between_obs <- projection_matrix(Z, tol = tol)$H
  H_within <- diag(n) - H_between_obs

  cls <- .classify_split_plot_terms(design, group, split$covariates, H_within, tol)
  if (length(cls$partial_terms)) {
    stop("Fixed term(s) ", paste(cls$partial_terms, collapse = ", "),
         " are neither constant nor fully crossed within '", group_label,
         "' (partial confounding with the random unit); the split-plot ",
         "denominator is not identified.", call. = FALSE)
  }
  X <- design$X
  assign <- cls$assign
  labels <- cls$labels
  intercept <- which(assign == 0L)
  between_terms <- cls$between_terms
  within_terms <- cls$within_terms
  between_columns <- cls$between_columns
  within_columns <- cls$within_columns

  if (!is.null(seed)) set.seed(seed)
  tested <- list()

  # ---- Between-unit (whole-plot) tests on equal-weight unit centroids -------
  if (length(between_terms)) {
    agg <- .aggregate_gower_units(design, group, tol = tol)
    Gc <- agg$G
    Xc <- X[agg$representatives, , drop = FALSE]        # unit-level design rows
    size_strata <- factor(agg$sizes, levels = sort(unique(agg$sizes)))
    plan_between <- .stratified_permutation_plan(size_strata, permutations, max_exact)

    stratum_reduced <- projection_matrix(
      Xc[, c(intercept, between_columns), drop = FALSE], tol = tol
    )
    stratum_H <- diag(n_units) - stratum_reduced$H
    stratum_df <- n_units - stratum_reduced$rank
    if (stratum_df <= 0L) {
      stop("No between-unit residual degrees of freedom remain for the ",
           "whole-plot denominator.", call. = FALSE)
    }
    stratum_ss <- trace_hgh(stratum_H, Gc)
    stratum_ms <- stratum_ss / stratum_df

    for (index in seq_along(between_terms)) {
      term <- between_terms[[index]]
      adjust <- between_terms[.ss_type_adjust(index, between_terms, hypothesis)]
      adjust_cols <- which(assign %in% match(adjust, labels))
      term_cols <- which(assign == match(term, labels))
      reduced <- projection_matrix(Xc[, c(intercept, adjust_cols), drop = FALSE], tol = tol)
      full <- projection_matrix(Xc[, c(intercept, adjust_cols, term_cols), drop = FALSE], tol = tol)
      num_H <- full$H - reduced$H
      num_df <- full$rank - reduced$rank
      if (num_df <= 0L) {
        stop("Whole-plot term `", term, "` is not estimable for Type ",
             hypothesis, " inference.", call. = FALSE)
      }
      term_basis <- (diag(n_units) - reduced$H) %*% Xc[, term_cols, drop = FALSE]
      if (.stratified_treatment_rank(term_basis, size_strata, tol) < num_df) {
        stop("Whole-plot term `", term, "` is not fully spanned within exact ",
             "unit-size permutation strata; its conditional permutation null ",
             "is not identifiable.", call. = FALSE)
      }
      fl <- .unit_freedman_lane(Gc, num_H, num_df, reduced$H, stratum_H,
                                stratum_df, plan_between, tol,
                                label = paste0("Permuting whole units for term `",
                                               term, "`"))
      tested[[term]] <- list(
        term = term, kind = "whole_plot", level = 1L,
        numerator_ss = fl$numerator_ss, numerator_df = num_df,
        numerator_ms = fl$numerator_ms, numerator_H = num_H,
        denominator_ms = stratum_ms, denominator_df = stratum_df,
        denominator_name = paste0(group_label, " stratum"),
        statistic = fl$observed, p_value = fl$p_value,
        exchangeable = "whole units (size-restricted)",
        permutations_used = nrow(plan_between$orders),
        attainable_permutations = plan_between$attainable,
        exact_permutation = plan_between$exact && nrow(plan_between$orders) > 0L,
        permutation_statistics = if (keep_permutations) fl$statistics else NULL
      )
    }
  } else {
    stratum_ss <- stratum_df <- stratum_ms <- NULL
    agg <- NULL
  }

  # ---- Within-unit (sub-plot) tests over the within-unit residual -----------
  within_proj <- .projection_in_subspace(H_within, X[, within_columns, drop = FALSE], tol = tol)
  residual_H <- H_within - within_proj$H
  residual_df <- (n - n_units) - within_proj$rank
  if (residual_df <= 0L) {
    stop("No within-unit residual degrees of freedom remain after the ",
         "within-unit terms.", call. = FALSE)
  }
  residual_ss <- trace_hgh(residual_H, G)
  residual_ms <- residual_ss / residual_df
  unit_blocks <- unname(split(seq_len(n), group))

  for (index in seq_along(within_terms)) {
    term <- within_terms[[index]]
    adjust <- within_terms[.ss_type_adjust(index, within_terms, hypothesis)]
    adjust_cols <- which(assign %in% match(adjust, labels))
    term_cols <- which(assign == match(term, labels))
    within_adjust <- if (length(adjust_cols)) H_within %*% X[, adjust_cols, drop = FALSE] else matrix(0, n, 0L)
    within_term <- H_within %*% X[, term_cols, drop = FALSE]
    reduced <- projection_matrix(within_adjust, tol = tol)
    full <- projection_matrix(cbind(within_adjust, within_term), tol = tol)
    num_H <- full$H - reduced$H
    num_df <- full$rank - reduced$rank
    if (num_df <= 0L) {
      stop("Within-unit term `", term, "` is not estimable for Type ",
           hypothesis, " inference.", call. = FALSE)
    }
    num_ss <- trace_hgh(num_H, G)
    num_ms <- num_ss / num_df
    reduced_full_H <- H_between_obs + reduced$H
    fl <- .within_unit_freedman_lane(G, num_H, num_df, reduced_full_H, residual_H,
                                     residual_df, num_ms / residual_ms, group,
                                     permutations, max_exact, tol)
    tested[[term]] <- list(
      term = term, kind = "within_unit", level = 2L,
      numerator_ss = num_ss, numerator_df = num_df, numerator_ms = num_ms,
      numerator_H = num_H, denominator_ms = residual_ms, denominator_df = residual_df,
      denominator_name = "Residual", statistic = num_ms / residual_ms,
      p_value = fl$p_value, exchangeable = "observations within units",
      permutations_used = fl$used,
      attainable_permutations = fl$attainable,
      exact_permutation = fl$exact,
      permutation_statistics = if (keep_permutations) fl$statistics else NULL
    )
  }

  specs <- stats::setNames(lapply(tested, function(item) {
    spec <- list(
      term = item$term, kind = item$kind,
      numerator = list(ss = item$numerator_ss, df = item$numerator_df,
                       ms = item$numerator_ms, projector = item$numerator_H),
      denominator = list(expression = paste0("MS[", item$denominator_name, "]"),
                         ms = item$denominator_ms, df = item$denominator_df),
      statistic = item$statistic, p_value = item$p_value,
      exchangeable_units = item$exchangeable,
      permutations_used = item$permutations_used,
      attainable_permutations = item$attainable_permutations,
      exact_permutation = item$exact_permutation,
      permutation_statistics = item$permutation_statistics
    )
    class(spec) <- "permanova_test_spec"
    spec
  }), names(tested))

  table <- do.call(rbind, lapply(specs, function(spec) data.frame(
    term = spec$term, df = spec$numerator$df, ss = spec$numerator$ss,
    ms = spec$numerator$ms, denominator = spec$denominator$expression,
    denominator_df = spec$denominator$df, pseudo_f = spec$statistic,
    p_perm = spec$p_value, stringsAsFactors = FALSE
  )))
  if (length(between_terms)) {
    table <- rbind(table, data.frame(
      term = paste0(group_label, " stratum"), df = stratum_df, ss = stratum_ss,
      ms = stratum_ms, denominator = NA_character_, denominator_df = NA_real_,
      pseudo_f = NA_real_, p_perm = NA_real_))
  }
  table <- rbind(table, data.frame(
    term = "Residual", df = residual_df, ss = residual_ss, ms = residual_ms,
    denominator = NA_character_, denominator_df = NA_real_, pseudo_f = NA_real_,
    p_perm = NA_real_))
  rownames(table) <- NULL

  out <- list(
    call = match.call(), design = design, table = table, tests = specs,
    hierarchy = group_label, covariates = split$covariates, factors = split$factors,
    supported_design = paste(
      "unequally replicated single random intercept with multiple categorical",
      "factors; whole-plot effects on equal-weight unit centroids, sub-plot",
      "effects over the within-unit residual"
    ),
    permutation = list(
      algorithm = "size-restricted whole-unit and within-unit residual permutation",
      used = max(vapply(specs, `[[`, numeric(1), "permutations_used")),
      exact = all(vapply(specs, `[[`, logical(1), "exact_permutation")), seed = seed
    )
  )
  class(out) <- "permanova_nested"
  out
}

# Freedman-Lane over observations within units: permute the reduced-model
# residual only inside each unit block (sub-plot exchangeability).
.within_unit_freedman_lane <- function(G, H_num, num_df, reduced_H, H_resid,
                                       resid_df, observed, unit,
                                       permutations, max_exact, tol,
                                       label = NULL) {
  plan <- .stratified_permutation_plan(unit, permutations, max_exact)
  stats_perm <- .freedman_lane_statistics(
    G, reduced_H, H_num, num_df, H_resid, resid_df, plan$orders,
    label = label, tol = tol
  )
  p_value <- NA_real_
  if (length(stats_perm)) {
    extreme <- sum(stats_perm >= observed - tol, na.rm = TRUE)
    p_value <- if (plan$exact) extreme / length(stats_perm) else
      (extreme + 1) / (length(stats_perm) + 1)
  }
  list(p_value = p_value, statistics = stats_perm, used = nrow(plan$orders),
       attainable = plan$attainable, exact = plan$exact && nrow(plan$orders) > 0L)
}
