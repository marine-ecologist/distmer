# Reduced-model term indices for a Type I/II/III hypothesis, computed within a
# set of fixed-term labels (a between-unit or within-unit block). Type I adjusts
# for the preceding terms; Type II for every term that does not *contain* the
# target (marginality); Type III for every other term. Returns positions in
# `term_labels`.
.ss_type_adjust <- function(target_position, term_labels, type) {
  if (identical(type, "I")) return(seq_len(target_position - 1L))
  if (identical(type, "III")) return(setdiff(seq_along(term_labels), target_position))
  target_vars <- strsplit(term_labels[[target_position]], ":", fixed = TRUE)[[1L]]
  containing <- vapply(seq_along(term_labels), function(j) {
    if (j == target_position) return(FALSE)
    other <- strsplit(term_labels[[j]], ":", fixed = TRUE)[[1L]]
    all(target_vars %in% other) && length(other) > length(target_vars)
  }, logical(1))
  setdiff(seq_along(term_labels), c(target_position, which(containing)))
}

#' Covariate-adjusted single-random-intercept PERMANOVA (split-plot / RCB ANCOVA)
#'
#' Tests categorical factors and one or more continuous covariates in a balanced
#' `(1 | group)` random-intercept design. Each fixed term is classified by
#' whether it varies within the grouping unit. Between-unit (whole-plot) factors,
#' constant within each unit, are tested over the between-unit random-intercept
#' stratum, exactly as in [permanova_mixed()]. Within-unit terms -- sub-plot or
#' randomized-block factors, observation-level covariates, and covariate-by-
#' factor slope interactions -- are partitioned in the within-unit space and
#' tested over the within-unit residual. This covers the distance-based
#' split-plot and randomized-block analysis of covariance.
#'
#' A numeric term that is constant within groups is assigned to the between-unit
#' stratum. A covariate that varies within groups is tested there while its
#' between-unit component is retained as nuisance adjustment for whole-plot
#' terms. [decompose_covariate()] may be used when separate between- and
#' within-unit slopes are required.
#'
#' @param design Canonical Design Object with one `(1 | group)` random term and
#'   a fixed model containing at least one continuous covariate.
#' @param permutations A non-negative integer or `"all"`.
#' @param seed Optional random seed.
#' @param keep_permutations Retain permutation statistics on each test spec.
#' @param max_exact Maximum permutation space allowed for exact enumeration.
#' @param tol Numerical tolerance.
#' @param hypothesis Sum-of-squares partition for the fixed terms: `"I"`
#'   (sequential), `"II"` (marginal, respecting containment), or `"III"`
#'   (partial). Types coincide for an orthogonal balanced fixed structure and
#'   differ only when fixed terms are non-orthogonal (interactions). Type III
#'   with interactions requires sum-to-zero contrasts, which the `permanova()`
#'   front door sets automatically.
#' @return An object of class `permanova_nested` (shared table/summary shape).
#' @export
permanova_mixed_covariate <- function(design, permutations = 999L, seed = NULL,
                                      hypothesis = "I", keep_permutations = FALSE,
                                      max_exact = 100000L, tol = 1e-10) {
  hypothesis <- match.arg(hypothesis, c("I", "II", "III"))
  if (!inherits(design, "permanova_design") || is.null(design$G)) {
    stop("`design` must be a compiled Design Object containing `distance`.",
         call. = FALSE)
  }
  if (length(design$random_terms) != 1L ||
      design$random_terms[[1L]]$origin != "main") {
    stop("Covariate-adjusted mixed inference requires one `(1 | group)` term.",
         call. = FALSE)
  }
  if (identical(hypothesis, "III")) {
    .validate_type3_contrasts(design, tol = tol)
  }
  split <- .nested_covariate_split(design)
  if (!length(split$covariates)) {
    stop("No continuous covariate found; use `permanova_mixed()`.", call. = FALSE)
  }

  group_label <- design$random_terms[[1L]]$label
  group <- droplevels(factor(design$groups[[group_label]]))
  n <- design$n
  G <- design$G
  n_units <- nlevels(group)
  unit_sizes <- tabulate(group)
  if (length(unique(unit_sizes)) != 1L) {
    stop("Balanced covariate-adjusted mixed inference requires equal ",
         "observations per group.", call. = FALSE)
  }
  if (unit_sizes[[1L]] < 2L) {
    stop("At least two observations per group are required.", call. = FALSE)
  }

  # Between-unit (group-mean) space and its within-unit complement.
  Z <- stats::model.matrix(~ 0 + group)
  between <- projection_matrix(Z, tol = tol)
  H_between <- between$H
  H_within <- diag(n) - H_between

  X <- design$X
  assign <- attr(X, "assign")
  labels <- attr(stats::terms(design$fixed_formula), "term.labels")
  intercept <- which(assign == 0L)
  term_variables <- lapply(labels, function(term) strsplit(term, ":", fixed = TRUE)[[1L]])
  term_has_covariate <- vapply(term_variables, function(v) any(v %in% split$covariates),
                               logical(1))
  # Classify each fixed term as between-unit (constant within every grouping
  # unit -> whole-plot, tested over the group stratum) or within-unit (varies
  # within units -> sub-plot factor, covariate, or slope interaction, tested
  # over the within-unit residual). Any covariate-involving term is within-unit.
  within_term_projection <- lapply(seq_along(labels), function(j) {
    cols <- which(assign == j)
    .projection_in_subspace(H_within, X[, cols, drop = FALSE], tol = tol)
  })
  varies_within_term <- vapply(within_term_projection, function(x) x$rank > 0L,
                               logical(1))
  is_within <- varies_within_term
  between_terms <- labels[!is_within]
  within_terms <- labels[is_within]
  between_columns <- which(assign %in% which(!is_within))
  within_columns <- which(assign %in% which(is_within))
  spanning_covariate_columns <- which(
    assign %in% which(term_has_covariate & is_within)
  )

  # Between-unit random-intercept stratum: the denominator for whole-plot
  # factors. The covariate's between-unit component is removed here too, so a
  # covariate with between-unit structure does not pollute the factor
  # denominator (a split-plot ANCOVA carries the covariate in both strata).
  between_covariate <- .projection_in_subspace(
    H_between, X[, spanning_covariate_columns, drop = FALSE], tol = tol
  )$projected
  between_fixed <- projection_matrix(
    cbind(X[, c(intercept, between_columns), drop = FALSE], between_covariate),
    tol = tol
  )
  stratum_H <- H_between - between_fixed$H
  stratum_df <- n_units - between_fixed$rank
  if (stratum_df <= 0L) {
    stop("No between-unit residual degrees of freedom remain for the factor ",
         "denominator.", call. = FALSE)
  }
  stratum_ss <- trace_hgh(stratum_H, G)
  stratum_ms <- stratum_ss / stratum_df

  # Within-unit residual after removing every within-unit term's component.
  within_projected <- .projection_in_subspace(
    H_within, X[, within_columns, drop = FALSE], tol = tol
  )
  residual_H <- H_within - within_projected$H
  residual_df <- (n - n_units) - within_projected$rank
  if (residual_df <= 0L) {
    stop("No within-unit residual degrees of freedom remain after the ",
         "within-unit terms.", call. = FALSE)
  }
  residual_ss <- trace_hgh(residual_H, G)
  residual_ms <- residual_ss / residual_df

  if (!is.null(seed)) set.seed(seed)
  tol_zero <- tol

  tested <- list()
  # Between-unit (whole-plot) factor tests are carried out entirely in the
  # between-unit space, adjusting for the covariate's *between-unit* component
  # (its own between-unit slope, which differs from the pooled slope) and any
  # preceding between-unit factors. The denominator is the group stratum.
  for (index in seq_along(between_terms)) {
    term <- between_terms[[index]]
    term_index <- match(term, labels)
    adjust_names <- between_terms[.ss_type_adjust(index, between_terms, hypothesis)]
    adjust_cols <- which(assign %in% match(adjust_names, labels))
    term_cols <- which(assign == term_index)
    reduced <- projection_matrix(
      cbind(X[, c(intercept, adjust_cols), drop = FALSE], between_covariate),
      tol = tol
    )
    full <- projection_matrix(
      cbind(X[, c(intercept, adjust_cols, term_cols), drop = FALSE],
            between_covariate),
      tol = tol
    )
    numerator_H <- full$H - reduced$H
    numerator_df <- full$rank - reduced$rank
    if (numerator_df <= 0L) {
      stop("Between-unit term `", term, "` is not estimable for Type ",
           hypothesis, " inference.", call. = FALSE)
    }
    numerator_ss <- trace_hgh(numerator_H, G)
    numerator_ms <- numerator_ss / numerator_df
    tested[[term]] <- list(
      term = term,
      kind = if (term_has_covariate[[term_index]]) "between-covariate" else "factor",
      level = 1L, reduced_H = reduced$H,
      numerator_H = numerator_H, numerator_df = numerator_df,
      numerator_ss = numerator_ss, numerator_ms = numerator_ms,
      denominator_H = stratum_H, denominator_df = stratum_df,
      denominator_ss = stratum_ss, denominator_ms = stratum_ms,
      denominator_name = "between-group stratum",
      statistic = numerator_ms / stratum_ms
    )
  }
  # Within-unit tests (sub-plot factors, covariates, and covariate-by-factor
  # slope interactions) are partitioned sequentially in the within-unit space
  # and tested over the within-unit residual.
  deepest_group_H <- H_between
  for (index in seq_along(within_terms)) {
    term <- within_terms[[index]]
    term_cols <- which(assign == match(term, labels))
    adjust_names <- within_terms[.ss_type_adjust(index, within_terms, hypothesis)]
    adjust_cols <- which(assign %in% match(adjust_names, labels))
    within_term <- .projection_in_subspace(
      H_within, X[, term_cols, drop = FALSE], tol = tol
    )$projected
    within_adjust <- if (length(adjust_cols)) {
      .projection_in_subspace(
        H_within, X[, adjust_cols, drop = FALSE], tol = tol
      )$projected
    } else matrix(0, n, 0L)
    full <- projection_matrix(cbind(within_adjust, within_term), tol = tol)
    reduced <- projection_matrix(within_adjust, tol = tol)
    numerator_H <- full$H - reduced$H
    numerator_df <- full$rank - reduced$rank
    if (numerator_df <= 0L) {
      stop("Within-unit term `", term, "` is not estimable for Type ",
           hypothesis, " inference.", call. = FALSE)
    }
    numerator_ss <- trace_hgh(numerator_H, G)
    numerator_ms <- numerator_ss / numerator_df
    kind <- if (term %in% labels[term_has_covariate]) "covariate" else "within-factor"
    tested[[term]] <- list(
      term = term, kind = kind, level = 2L,
      reduced_H = deepest_group_H + reduced$H,
      numerator_H = numerator_H, numerator_df = numerator_df,
      numerator_ss = numerator_ss, numerator_ms = numerator_ms,
      denominator_H = residual_H, denominator_df = residual_df,
      denominator_ss = residual_ss, denominator_ms = residual_ms,
      denominator_name = "Residual",
      statistic = numerator_ms / residual_ms
    )
  }

  specs <- .mixed_covariate_permutation(
    tested, design, group, H_within, permutations, seed, max_exact,
    keep_permutations, tol_zero
  )

  table <- do.call(rbind, lapply(specs, function(spec) data.frame(
    term = spec$term, df = spec$numerator$df, ss = spec$numerator$ss,
    ms = spec$numerator$ms, denominator = spec$denominator$expression,
    denominator_df = spec$denominator$df, pseudo_f = spec$statistic,
    p_perm = spec$p_value, stringsAsFactors = FALSE
  )))
  table <- rbind(
    table,
    data.frame(term = paste0(group_label, " (between-group stratum)"),
               df = stratum_df, ss = stratum_ss, ms = stratum_ms,
               denominator = NA_character_, denominator_df = NA_real_,
               pseudo_f = NA_real_, p_perm = NA_real_),
    data.frame(term = "Residual", df = residual_df, ss = residual_ss,
               ms = residual_ms, denominator = NA_character_,
               denominator_df = NA_real_, pseudo_f = NA_real_, p_perm = NA_real_)
  )
  rownames(table) <- NULL

  # Variance components: between-unit and within-unit trace components.
  m <- unit_sizes[[1L]]
  between_component <- (stratum_ms - residual_ms) / m
  variance_components <- list(
    unrestricted = c(residual = residual_ms, `group` = between_component),
    method = "balanced split-plot trace method of moments (covariate-adjusted residual)"
  )

  out <- list(
    call = match.call(), design = design, table = table, tests = specs,
    variance_components = variance_components,
    component_labels = stats::setNames(group_label, group_label),
    hierarchy = group_label, covariates = split$covariates,
    factors = split$factors,
    supported_design = paste(
      "balanced single random intercept with categorical factors adjusted for",
      "an observation-level continuous covariate (split-plot ANCOVA)"
    ),
    permutation = list(
      algorithm = "term-specific Freedman-Lane permutation, covariate held in the reduced model",
      used = max(vapply(specs, `[[`, numeric(1), "permutations_used")),
      seed = seed
    )
  )
  class(out) <- "permanova_nested"
  out
}

# Level-specific Freedman-Lane permutation for the split-plot ANCOVA: factors
# permute whole groups (between-unit), covariate terms permute observations
# within groups (within-unit), each holding its reduced model fixed.
.mixed_covariate_permutation <- function(tested, design, group, H_within,
                                         permutations, seed, max_exact,
                                         keep_permutations, tol) {
  n <- design$n
  group_index <- as.integer(group)
  unit_blocks <- split(seq_len(n), group)
  stats::setNames(lapply(tested, function(item) {
    resid <- diag(n) - item$reduced_H
    plan <- .mixed_covariate_plan(item$kind, group, permutations, max_exact)
    permuted <- numeric(nrow(plan$orders))
    if (nrow(plan$orders)) {
      for (b in seq_len(nrow(plan$orders))) {
        idx <- plan$index(plan$orders[b, ])
        A <- item$reduced_H + resid[idx, , drop = FALSE]
        G_perm <- A %*% design$G %*% t(A)
        num <- trace_hgh(item$numerator_H, G_perm) / item$numerator_df
        den <- trace_hgh(item$denominator_H, G_perm) / item$denominator_df
        permuted[[b]] <- num / den
      }
    }
    p_value <- NA_real_
    if (length(permuted)) {
      extreme <- sum(permuted >= item$statistic - tol, na.rm = TRUE)
      p_value <- if (plan$exact) extreme / length(permuted) else
        (extreme + 1) / (length(permuted) + 1)
    }
    spec <- list(
      term = item$term, kind = item$kind,
      numerator = list(ss = item$numerator_ss, df = item$numerator_df,
                       ms = item$numerator_ms, projector = item$numerator_H),
      denominator = list(
        expression = paste0("MS[", item$denominator_name, "]"),
        ss = item$denominator_ss, df = item$denominator_df,
        ms = item$denominator_ms, projector = item$denominator_H),
      statistic = item$statistic, p_value = p_value,
      exchangeable_units = plan$level_name,
      permutation_algorithm = plan$algorithm,
      attainable_permutations = plan$attainable,
      permutations_used = nrow(plan$orders),
      exact_permutation = plan$exact && nrow(plan$orders) > 0L,
      permutation_statistics = if (keep_permutations) permuted else NULL
    )
    class(spec) <- "permanova_test_spec"
    spec
  }), names(tested))
}

# Whole-group permutation index for factors, within-group for covariates.
.mixed_covariate_plan <- function(kind, group, permutations, max_exact) {
  n <- length(group)
  units <- unname(split(seq_len(n), group))
  n_units <- length(units)
  if (kind %in% c("factor", "between-covariate")) {
    unit_plan <- .permutation_orders(n_units, permutations, max_exact)
    orders <- if (nrow(unit_plan$orders)) {
      t(vapply(seq_len(nrow(unit_plan$orders)), function(i) {
        .block_permutation_index(units, unit_plan$orders[i, ])
      }, integer(n)))
    } else matrix(integer(), 0L, n)
    list(orders = orders, index = identity, exact = unit_plan$exact,
         attainable = unit_plan$attainable, level_name = "whole groups",
         algorithm = "whole-group residual permutation")
  } else {
    within_plan <- .stratified_permutation_plan(
      group, permutations, max_exact
    )
    list(orders = within_plan$orders, index = identity,
         exact = within_plan$exact, attainable = within_plan$attainable,
         level_name = "within groups",
         algorithm = "within-group residual permutation")
  }
}
