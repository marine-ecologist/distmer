#' Split a covariate into between-unit and within-unit components
#'
#' Separates one or more continuous covariates into a between-unit part (the
#' grouping-unit mean) and a within-unit part (the deviation from that mean).
#' Entering both parts as predictors lets engines that support both strata
#' adjust the between-unit and within-unit relationships separately. The within
#' part is a genuine observation-level covariate suitable for
#' [permanova_nested_covariate()]; the between part is a unit-level summary and
#' must only be entered in an engine that assigns numeric terms to that level
#' (for example [permanova_mixed_covariate()]).
#'
#' @param data A data frame.
#' @param covariate Character vector of continuous covariate column names.
#' @param by Character vector of grouping column name(s) defining the units; the
#'   between/within split is taken with respect to their combination.
#' @param center Grand-centre the between-unit component (default `TRUE`) so it
#'   is a deviation from the overall mean. The within-unit component always has
#'   mean zero within every unit regardless of this setting.
#' @param suffix Length-two character vector of column-name suffixes for the
#'   between and within parts.
#' @return `data` with two added columns per covariate. With `center = FALSE`
#'   the between and within columns sum to the original covariate.
#' @examples
#' df <- data.frame(site = factor(rep(1:3, each = 4)), x = rnorm(12))
#' decompose_covariate(df, "x", by = "site")
#' @export
decompose_covariate <- function(data, covariate, by, center = TRUE,
                                suffix = c("_between", "_within")) {
  data <- as.data.frame(data)
  if (!is.character(covariate) || !length(covariate) ||
      length(setdiff(covariate, names(data)))) {
    stop("`covariate` must name numeric columns of `data`.", call. = FALSE)
  }
  if (!is.character(by) || !length(by) || length(setdiff(by, names(data)))) {
    stop("`by` must name one or more grouping columns of `data`.", call. = FALSE)
  }
  if (length(suffix) != 2L || !is.character(suffix) || anyNA(suffix) ||
      identical(suffix[[1L]], suffix[[2L]])) {
    stop("`suffix` must be two distinct strings.", call. = FALSE)
  }
  groups <- interaction(lapply(data[by], factor), drop = TRUE, lex.order = TRUE)
  for (name in covariate) {
    value <- data[[name]]
    if (!is.numeric(value) || is.factor(value)) {
      stop("Covariate `", name, "` must be an ordinary numeric column.",
           call. = FALSE)
    }
    unit_mean <- stats::ave(value, groups, FUN = function(z) mean(z))
    within <- value - unit_mean
    between <- if (isTRUE(center)) unit_mean - mean(value) else unit_mean
    data[[paste0(name, suffix[[1L]])]] <- between
    data[[paste0(name, suffix[[2L]])]] <- within
  }
  data
}

# PRIMER-style fixed-effect "components" for a covariate-adjusted nested fit:
# for each fixed term, (MS_term - MS_denominator) / c, where c = tr(H Z Z') / df.
# These are effect sizes (they can be negative), NOT random-effect variances, so
# summary() prints them in a separate, clearly labelled block.
.nested_covariate_fixed_components <- function(object) {
  fit <- object$fit
  design <- object$design
  if (is.null(fit$tests)) return(NULL)
  labels <- attr(stats::terms(design$fixed_formula), "term.labels")
  assignment <- attr(design$X, "assign")
  rows <- lapply(fit$tests, function(spec) {
    term <- spec$term
    variables <- strsplit(term, ":", fixed = TRUE)[[1L]]
    if (identical(spec$kind, "factor")) {
      cells <- interaction(lapply(design$data[variables], factor), drop = TRUE)
      basis <- stats::model.matrix(~ 0 + cells)
    } else {
      basis <- design$X[, assignment == match(term, labels), drop = FALSE]
    }
    coefficient <- trace_product(spec$numerator$projector, tcrossprod(basis)) /
      spec$numerator$df
    estimate <- if (abs(coefficient) > 1e-10) {
      (spec$numerator$ms - spec$denominator$ms) / coefficient
    } else NA_real_
    data.frame(
      Source = paste0("S(", term, ")"), Estimate = estimate,
      Sq.root = sign(estimate) * sqrt(abs(estimate)),
      check.names = FALSE, stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

.nested_covariate_split <- function(design) {
  fixed_variables <- all.vars(design$fixed_formula)
  is_covariate <- vapply(fixed_variables, function(variable) {
    value <- design$data[[variable]]
    is.numeric(value) && !is.factor(value)
  }, logical(1))
  list(
    covariates = fixed_variables[is_covariate],
    factors = fixed_variables[!is_covariate]
  )
}

#' Covariate-adjusted balanced nested PERMANOVA
#'
#' Tests categorical factors in a balanced `/`-generated nested hierarchy while
#' adjusting for one or more continuous covariates. Each factor is tested for an
#' effect *beyond* the covariate relationship: its numerator is the sum of
#' squares added over the covariate(s), and its denominator is the same exact
#' nested expected-mean-square stratum used by [permanova_nested()]. The
#' covariate is treated as an observation-level nuisance (charged to the
#' residual) and is itself tested over that residual.
#'
#' A single-slope covariate cannot be split cleanly across nested strata, so a
#' covariate with strong between-unit structure should be decomposed by the user
#' into explicit between- and within-unit predictors before testing.
#'
#' @param design Canonical Design Object whose fixed model contains at least one
#'   categorical factor and at least one continuous covariate, with one balanced
#'   `/`-generated random hierarchy.
#' @param permutations A non-negative integer or `"all"`.
#' @param seed Optional random seed.
#' @param hypothesis Sum-of-squares partition for the fixed factors: `"I"`
#'   (sequential), `"II"` (marginal), or `"III"` (partial). Types coincide for an
#'   orthogonal balanced fixed structure; Type III with interactions requires
#'   sum-to-zero contrasts, set automatically by the `permanova()` front door.
#' @param keep_permutations Retain permutation statistics.
#' @param max_exact Maximum exact permutation space.
#' @param tol Numerical tolerance.
#' @return An object of class `permanova_nested`.
#' @export
permanova_nested_covariate <- function(design, permutations = 999L, seed = NULL,
                                       hypothesis = "I", keep_permutations = FALSE,
                                       max_exact = 100000L, tol = 1e-10) {
  hypothesis <- match.arg(hypothesis, c("I", "II", "III"))
  if (identical(hypothesis, "III")) {
    .validate_type3_contrasts(design, tol = tol)
  }
  split <- .nested_covariate_split(design)
  if (!length(split$covariates)) {
    stop("No continuous covariate found; use `permanova_nested()`.", call. = FALSE)
  }
  if (!length(split$factors)) {
    stop("At least one categorical factor is required.", call. = FALSE)
  }

  # Factor-only nested structure: reuse the exact-EMS validator and its strata.
  nesting <- Filter(function(x) x$origin == "nesting", design$random_terms)
  if (!length(nesting)) {
    stop("Covariate-adjusted nested inference requires a `/`-generated hierarchy.",
         call. = FALSE)
  }
  chain <- nesting[[length(nesting)]]$expression
  factor_formula <- stats::as.formula(
    paste("~", paste(c(split$factors, sprintf("(1 | %s)", chain)), collapse = " + ")),
    env = environment(design$formula)
  )
  factor_design <- compile_permanova_design(
    factor_formula, design$data, design$distance, tol = tol
  )
  structure <- .validate_nested_design(factor_design, tol = tol)
  if (!is.null(seed)) set.seed(seed)

  X <- design$X
  assign <- attr(X, "assign")
  labels <- attr(stats::terms(design$fixed_formula), "term.labels")
  intercept <- which(assign == 0L)
  # A term is covariate-involving if any of its variables is a covariate; this
  # captures the covariate main effect and any covariate-by-factor slope
  # interaction (the homogeneity-of-slopes terms). Pure-factor terms are tested
  # over their nested strata; covariate-involving terms are observation-level
  # and tested over the residual.
  term_variables <- lapply(labels, function(term) {
    strsplit(term, ":", fixed = TRUE)[[1L]]
  })
  term_has_covariate <- vapply(term_variables, function(variables) {
    any(variables %in% split$covariates)
  }, logical(1))
  factor_terms <- labels[!term_has_covariate]
  covariate_terms <- labels[term_has_covariate]
  covariate_columns <- which(assign %in% which(term_has_covariate))

  # The covariate is an observation-level nuisance: it is charged to the deepest
  # within-unit stratum. Its within-unit component (deepest_H %*% C) is removed
  # from that stratum to form the residual, so degrees of freedom balance.
  deepest <- structure$strata[[length(structure$strata)]]
  covariate_matrix <- X[, covariate_columns, drop = FALSE]
  covariate_within <- .projection_in_subspace(
    deepest$H, covariate_matrix, tol = tol
  )
  for (term in covariate_terms) {
    term_cols <- which(assign == match(term, labels))
    term_within <- .projection_in_subspace(
      deepest$H, X[, term_cols, drop = FALSE], tol = tol
    )
    if (term_within$rank == 0L) {
      stop(
        "Numeric term `", term, "` has no estimable variation in the deepest ",
        "nested stratum. `permanova_nested_covariate()` supports observation-level ",
        "covariates only; omit the between-unit component or use an engine that ",
        "assigns it an explicit higher-level denominator.", call. = FALSE
      )
    }
  }
  residual_H <- deepest$H - covariate_within$H
  residual_df <- deepest$df - covariate_within$rank
  if (residual_df <= 0L) {
    stop("No within-unit residual degrees of freedom remain after the covariate.",
         call. = FALSE)
  }
  residual_ss <- trace_hgh(residual_H, design$G)
  residual_ms <- residual_ss / residual_df
  deepest_group_H <- diag(design$n) - deepest$H
  observation_level <- length(structure$level_names)

  tested <- list()
  # --- factor tests: numerator adjusted for covariate + preceding factors ---
  for (index in seq_along(factor_terms)) {
    term <- factor_terms[[index]]
    term_index <- match(term, labels)
    adjust_names <- factor_terms[.ss_type_adjust(index, factor_terms, hypothesis)]
    reduced_cols <- unique(c(
      intercept, covariate_columns,
      which(assign %in% match(adjust_names, labels))
    ))
    term_cols <- which(assign == term_index)
    full_cols <- unique(c(reduced_cols, term_cols))
    reduced <- projection_matrix(X[, reduced_cols, drop = FALSE], tol = tol)
    tested_proj <- projection_matrix(X[, full_cols, drop = FALSE], tol = tol)
    numerator_H <- tested_proj$H - reduced$H
    numerator_df <- tested_proj$rank - reduced$rank
    if (numerator_df <= 0L) {
      stop("Fixed term `", term, "` is not estimable after covariate adjustment.",
           call. = FALSE)
    }
    level <- structure$term_levels[[term]]
    stratum <- structure$strata[[level]]
    numerator_ss <- trace_hgh(numerator_H, design$G)
    numerator_ms <- numerator_ss / numerator_df
    statistic <- numerator_ms / stratum$ms

    tested[[term]] <- list(
      term = term, level = level, kind = "factor",
      numerator_H = numerator_H, numerator_df = numerator_df,
      reduced_H = reduced$H,
      numerator_ss = numerator_ss, numerator_ms = numerator_ms,
      denominator_H = stratum$H, denominator_df = stratum$df,
      denominator_ss = stratum$ss, denominator_ms = stratum$ms,
      denominator_name = stratum$term, statistic = statistic
    )
  }
  # --- covariate-involving terms (main effect and any slope interactions) ---
  # Each is tested over the within residual using its within-unit component and
  # the requested Type I/II/III reduced model.
  residual_name <- paste0("Within ", tail(structure$labels, 1L), " after covariate")
  for (index in seq_along(covariate_terms)) {
    term <- covariate_terms[[index]]
    term_cols <- which(assign == match(term, labels))
    adjust_names <- covariate_terms[
      .ss_type_adjust(index, covariate_terms, hypothesis)
    ]
    adjust_cols <- which(assign %in% match(adjust_names, labels))
    within_term <- .projection_in_subspace(
      deepest$H, X[, term_cols, drop = FALSE], tol = tol
    )$projected
    within_adjust <- if (length(adjust_cols)) {
      .projection_in_subspace(
        deepest$H, X[, adjust_cols, drop = FALSE], tol = tol
      )$projected
    } else matrix(0, design$n, 0L)
    tested_proj <- projection_matrix(cbind(within_adjust, within_term), tol = tol)
    reduced_proj <- projection_matrix(within_adjust, tol = tol)
    numerator_H <- tested_proj$H - reduced_proj$H
    numerator_df <- tested_proj$rank - reduced_proj$rank
    if (numerator_df <= 0L) {
      stop("Covariate term `", term, "` is not estimable for Type ", hypothesis,
           " inference in the deepest nested stratum.", call. = FALSE)
    }
    numerator_ss <- trace_hgh(numerator_H, design$G)
    numerator_ms <- numerator_ss / numerator_df
    tested[[term]] <- list(
      term = term, level = observation_level, kind = "covariate",
      numerator_H = numerator_H, numerator_df = numerator_df,
      reduced_H = deepest_group_H + reduced_proj$H,
      numerator_ss = numerator_ss, numerator_ms = numerator_ms,
      denominator_H = residual_H, denominator_df = residual_df,
      denominator_ss = residual_ss, denominator_ms = residual_ms,
      denominator_name = residual_name, statistic = numerator_ms / residual_ms
    )
  }

  # PRIMER-style containment labels for the strata (Site(Location),
  # Area(Site(Location)), Residual), shared by the table and denominators.
  strata_terms <- vapply(structure$strata, `[[`, character(1), "term")
  primer_labels <- .primer_nested_labels(
    design$data, structure$labels, strata_terms, factor_terms, residual_name
  )
  strata_map <- primer_labels$strata
  relabel <- function(name) if (name %in% names(strata_map)) strata_map[[name]] else name

  # --- level-specific Freedman-Lane permutation, covariate held in reduced ---
  # Both factors and the observation-level covariate permute residuals within
  # their assigned hierarchy level, moving descendants intact.
  level_plans <- list()
  specs <- setNames(lapply(names(tested), function(key) {
    item <- tested[[key]]
    level_key <- as.character(item$level)
    if (is.null(level_plans[[level_key]])) {
      level_plans[[level_key]] <<- .nested_permutation_plan(
        structure, item$level, permutations, max_exact
      )
    }
    plan <- level_plans[[level_key]]
    blocks <- plan$unit_blocks
    permuted <- .freedman_lane_statistics(
      design$G, item$reduced_H,
      item$numerator_H, item$numerator_df,
      item$denominator_H, item$denominator_df,
      plan$orders,
      index = function(order) .block_permutation_index(blocks, order),
      label = paste0("Permuting term `", item$term, "`"), tol = tol
    )
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
      reduced_model = list(projector = item$reduced_H),
      denominator = list(expression = paste0("MS[", relabel(item$denominator_name), "]"),
                         ss = item$denominator_ss, df = item$denominator_df,
                         ms = item$denominator_ms, projector = item$denominator_H),
      statistic = item$statistic, p_value = p_value,
      exchangeable_units = plan$level_name,
      permutation_algorithm = plan$algorithm,
      attainable_permutations = plan$attainable,
      permutations_used = nrow(plan$orders),
      adjusted_for = if (item$kind == "factor") {
        c(split$covariates, factor_terms[factor_terms != item$term])
      } else factor_terms,
      permutation_statistics = if (keep_permutations) permuted else NULL
    )
    class(spec) <- "permanova_test_spec"
    spec
  }), names(tested))
  names(specs) <- vapply(tested, `[[`, character(1), "term")

  table <- do.call(rbind, lapply(specs, function(spec) data.frame(
    term = spec$term, df = spec$numerator$df, ss = spec$numerator$ss,
    ms = spec$numerator$ms, denominator = spec$denominator$expression,
    denominator_df = spec$denominator$df, pseudo_f = spec$statistic,
    p_perm = spec$p_value, stringsAsFactors = FALSE
  )))
  for (index in seq_along(structure$strata)) {
    stratum <- structure$strata[[index]]
    # the deepest stratum is now the covariate-adjusted within-unit residual
    if (index == length(structure$strata)) {
      table <- rbind(table, data.frame(
        term = relabel(stratum$term), df = residual_df,
        ss = residual_ss, ms = residual_ms, denominator = NA_character_,
        denominator_df = NA_real_, pseudo_f = NA_real_, p_perm = NA_real_
      ))
    } else {
      table <- rbind(table, data.frame(
        term = relabel(stratum$term), df = stratum$df, ss = stratum$ss,
        ms = stratum$ms, denominator = NA_character_, denominator_df = NA_real_,
        pseudo_f = NA_real_, p_perm = NA_real_
      ))
    }
  }
  rownames(table) <- NULL

  # Nested variance components from the trace-derived EMS system, using the
  # covariate-adjusted residual as the deepest mean square.
  selected <- structure$strata[rev(names(structure$strata))]
  coefficient_matrix <- do.call(rbind, lapply(selected, `[[`, "ems"))
  observed_ms <- vapply(selected, `[[`, numeric(1), "ms")
  observed_ms[[1L]] <- residual_ms
  components <- tryCatch(
    stats::setNames(as.vector(solve(coefficient_matrix, observed_ms)),
                    structure$components),
    error = function(e) NULL
  )
  variance_components <- if (is.null(components)) NULL else list(
    unrestricted = components,
    method = "balanced trace-derived nested EMS equations (covariate-adjusted residual)"
  )

  out <- list(
    call = match.call(), design = design, table = table, tests = specs,
    variance_components = variance_components,
    component_labels = primer_labels$components,
    hierarchy = structure$labels, covariates = split$covariates,
    factors = split$factors,
    supported_design = paste(
      "balanced nested random hierarchy with categorical factors adjusted for",
      "an observation-level continuous covariate"
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

# PRIMER-style nesting labels for the strata of a balanced nested design.
# Turns the internal stratum names ("Site within fixed model", "Site:Area
# within Site and fixed model", "Within Site:Area after ...") into the
# containment notation "Site(Location)", "Area(Site(Location))", "Residual".
# The parent chain is the fixed factor(s) that are constant within the top
# random unit, then each nested random level in turn. Returns a named character
# vector mapping every internal stratum name (and the covariate residual name)
# to its PRIMER label; unknown names map to themselves.
.primer_nested_labels <- function(data, hierarchy, strata_terms, factor_terms,
                                  residual_name = NULL) {
  own <- character(length(hierarchy))
  own[[1L]] <- hierarchy[[1L]]
  for (i in seq_along(hierarchy)[-1L]) {
    parent_vars <- strsplit(hierarchy[[i - 1L]], ":", fixed = TRUE)[[1L]]
    this_vars <- strsplit(hierarchy[[i]], ":", fixed = TRUE)[[1L]]
    own[[i]] <- paste(setdiff(this_vars, parent_vars), collapse = ":")
  }
  top_vars <- strsplit(hierarchy[[1L]], ":", fixed = TRUE)[[1L]]
  top_unit <- interaction(data[top_vars], drop = TRUE)
  grouping <- Filter(function(term) {
    vars <- strsplit(term, ":", fixed = TRUE)[[1L]]
    if (length(vars) != 1L || !vars %in% names(data)) return(FALSE)
    values <- data[[vars]]
    if (!is.factor(values) && !is.character(values)) return(FALSE)
    all(tapply(as.character(values), top_unit,
               function(v) length(unique(v)) == 1L))
  }, factor_terms)
  chain <- character(length(own))
  parent <- paste(unlist(grouping), collapse = ", ")
  for (i in seq_along(own)) {
    chain[[i]] <- if (nzchar(parent)) paste0(own[[i]], "(", parent, ")") else own[[i]]
    parent <- chain[[i]]
  }
  labels <- stats::setNames(strata_terms, strata_terms)
  for (i in seq_along(chain)) labels[[strata_terms[[i]]]] <- chain[[i]]
  labels[[strata_terms[[length(strata_terms)]]]] <- "Residual"
  if (!is.null(residual_name)) labels[[residual_name]] <- "Residual"
  # component labels are keyed by the hierarchy variance-component names so the
  # ECV V() rows can carry the same containment notation as the table.
  list(strata = labels, components = stats::setNames(chain, hierarchy))
}
