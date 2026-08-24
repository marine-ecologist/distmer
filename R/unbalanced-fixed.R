
# `%||%`, .warnf, and .stopf are defined once in R/utils.R.


.unbalanced_term_variables <- function(term) {
  variables <- strsplit(term, ":", fixed = TRUE)[[1L]]
  if (!length(variables) || any(!nzchar(variables))) {
    stop("Unbalanced hypotheses require ordinary variable and `:` term labels.",
         call. = FALSE)
  }
  variables
}

.validate_type3_contrasts <- function(design, tol) {
  if (!any(grepl(":", design$fixed_effects, fixed = TRUE))) {
    return(invisible(TRUE))
  }
  variables <- unique(unlist(lapply(
    design$fixed_effects, .unbalanced_term_variables
  )))
  factor_variables <- variables[vapply(
    design$data[variables], is.factor, logical(1)
  )]
  invalid <- factor_variables[vapply(factor_variables, function(variable) {
    contrast <- stats::contrasts(design$data[[variable]])
    is.null(contrast) || any(abs(colSums(contrast)) > sqrt(tol))
  }, logical(1))]
  if (length(invalid)) {
    stop(
      "Type III hypotheses with interactions require explicit sum-to-zero contrasts such as `contr.sum` for: ",
      paste(invalid, collapse = ", "), ".",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.validate_unbalanced_fixed_design <- function(design, hypothesis, tol) {
  if (!inherits(design, "permanova_design")) {
    stop("`design` must be produced by `compile_permanova_design()`.",
         call. = FALSE)
  }
  if (is.null(design$G)) {
    stop("The Design Object must include `distance`.", call. = FALSE)
  }
  if (length(design$random_terms)) {
    stop(
      "Phase 10 unbalanced inference supports fixed effects only; unbalanced mixed-model EMS denominators are not inferred implicitly.",
      call. = FALSE
    )
  }
  if (!length(design$fixed_effects)) {
    stop("At least one fixed-effect term is required.", call. = FALSE)
  }
  if (design$estimability$fixed_rank >= design$n) {
    stop("No full-model residual degrees of freedom remain.", call. = FALSE)
  }
  if (hypothesis == "III") {
    if (design$estimability$fixed_rank < ncol(design$X)) {
      stop(
        "Type III inference requires a full-rank fixed model; aliased columns from incomplete cells are unsupported.",
        call. = FALSE
      )
    }
    .validate_type3_contrasts(design, tol = tol)
  }
  invisible(TRUE)
}

.projection_for_columns <- function(X, columns, tol) {
  projection_matrix(X[, columns, drop = FALSE], tol = tol)
}

#' Derive estimable hypotheses for an unbalanced fixed model
#'
#' Constructs sequential Type I, marginal Type II, or full-model Type III
#' projection hypotheses. Type II is restricted to terms without containing
#' higher-order interactions; Type III requires a full-rank model and explicit
#' sum-to-zero contrasts when interactions are present.
#'
#' @param design Canonical fixed-effect Design Object.
#' @param hypothesis One of `"I"`, `"II"`, or `"III"`.
#' @param terms Optional fixed terms to derive.
#' @param tol Numerical tolerance.
#' @return A named list of hypothesis projectors and estimability metadata.
#' @export
derive_unbalanced_hypotheses <- function(
    design, hypothesis = c("I", "II", "III"), terms = NULL,
    tol = 1e-10) {
  hypothesis <- match.arg(hypothesis)
  .validate_unbalanced_fixed_design(design, hypothesis, tol = tol)
  labels <- design$fixed_effects
  if (is.null(terms)) terms <- labels
  if (!is.character(terms) || !length(terms) || anyDuplicated(terms) ||
      length(setdiff(terms, labels))) {
    stop("`terms` must be unique fixed-effect term labels in the compiled design.",
         call. = FALSE)
  }
  assignment <- attr(design$X, "assign")
  if (is.null(assignment)) {
    stop("The compiled fixed-model matrix lacks term assignments.",
         call. = FALSE)
  }
  intercept_columns <- which(assignment == 0L)

  all_specs <- setNames(lapply(seq_along(labels), function(index) {
    target <- labels[[index]]
    target_variables <- .unbalanced_term_variables(target)
    term_columns <- which(assignment == index)
    nominal_df <- length(term_columns)
    containing <- which(vapply(seq_along(labels), function(other_index) {
      if (other_index == index) return(FALSE)
      other_variables <- .unbalanced_term_variables(labels[[other_index]])
      all(target_variables %in% other_variables) &&
        length(other_variables) > length(target_variables)
    }, logical(1)))

    if (hypothesis == "I") {
      previous_terms <- seq_len(index - 1L)
      reduced_columns <- c(
        intercept_columns,
        which(assignment %in% previous_terms)
      )
      test_columns <- c(reduced_columns, term_columns)
      adjusted_for <- labels[previous_terms]
      excluded <- character()
    } else if (hypothesis == "II") {
      # Type II respects marginality: a term is adjusted for every other term
      # that does not contain it, and the higher-order interactions that do
      # contain it are held out of the reduced model. SS(A | 1, B) for A in
      # A + B + A:B. This is estimable in the presence of interactions and
      # matches the SAS/`car`/PRIMER "Type II (conditional)" hypothesis.
      adjust_terms <- setdiff(seq_along(labels), c(index, containing))
      reduced_columns <- c(
        intercept_columns, which(assignment %in% adjust_terms)
      )
      test_columns <- c(reduced_columns, term_columns)
      adjusted_for <- labels[adjust_terms]
      excluded <- labels[containing]
    } else {
      reduced_columns <- which(assignment != index)
      test_columns <- seq_len(ncol(design$X))
      adjusted_for <- labels[-index]
      excluded <- character()
    }

    reduced <- .projection_for_columns(
      design$X, unique(reduced_columns), tol = tol
    )
    tested <- .projection_for_columns(
      design$X, unique(test_columns), tol = tol
    )
    df <- tested$rank - reduced$rank
    H <- tested$H - reduced$H
    projector_error <- if (df > 0L) max(abs(H %*% H - H)) else Inf
    fully_estimable <- df == nominal_df && df > 0L &&
      projector_error <= max(100 * tol, sqrt(.Machine$double.eps))
    list(
      term = target, hypothesis_type = hypothesis,
      H_reduced = reduced$H, H_tested = tested$H, H = H,
      reduced_rank = reduced$rank, tested_rank = tested$rank,
      df = df, nominal_df = nominal_df,
      estimable = df > 0L, fully_estimable = fully_estimable,
      adjusted_for = adjusted_for,
      excluded_higher_order = excluded,
      containing_interactions = labels[containing],
      aliased_columns = design$estimability$aliased,
      projector_error = projector_error,
      reason = if (!df) {
        paste0("Term `", target, "` adds no estimable dimension.")
      } else if (!fully_estimable) {
        paste0(
          "Term `", target, "` has ", df, " estimable df of ",
          nominal_df, " nominal df."
        )
      } else "fully estimable"
    )
  }), labels)
  all_specs[terms]
}

.validate_unbalanced_hypotheses <- function(specifications) {
  invalid <- specifications[!vapply(
    specifications, `[[`, logical(1), "fully_estimable"
  )]
  if (length(invalid)) {
    details <- vapply(invalid, `[[`, character(1), "reason")
    stop(
      "Selected unbalanced hypothesis is not fully estimable: ",
      paste(details, collapse = " "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.unbalanced_test_specification <- function(
    term, specification, design, denominator, plan, tol) {
  numerator_ss <- trace_hgh(specification$H, design$G)
  numerator_ms <- numerator_ss / specification$df
  denominator_ss <- trace_hgh(denominator$H, design$G)
  denominator_ms <- denominator_ss / denominator$df
  if (.mean_square_is_zero(denominator_ms, design$G, tol)) {
    stop("The full-model residual mean square is numerically zero.",
         call. = FALSE)
  }
  observed <- numerator_ms / denominator_ms
  statistics <- .freedman_lane_statistics(
    design$G, specification$H_reduced,
    specification$H, specification$df,
    denominator$H, denominator$df,
    plan$orders,
    label = paste0("Permuting term `", term, "`"), tol = tol
  )
  inference <- .permutation_inference(statistics, observed, plan$exact, tol)
  p_value <- inference$p_value
  notes <- character()
  if (!isTRUE(design$spectral$euclidean)) {
    notes <- c(
      notes,
      "Empirical inference retained the complete non-Euclidean Gower spectrum."
    )
  }
  if (numerator_ss < -tol || denominator_ss < -tol) {
    notes <- c(notes, "A negative distance-based sum of squares was retained.")
  }
  out <- list(
    term = term,
    hypothesis = paste0(
      "Type ", specification$hypothesis_type,
      " estimable function for `", term, "` is zero"
    ),
    hypothesis_type = specification$hypothesis_type,
    numerator = list(
      ss = numerator_ss, df = specification$df, ms = numerator_ms,
      projector = specification$H
    ),
    denominator = list(
      expression = "MS[full-model Residual]",
      ss = denominator_ss, df = denominator$df, ms = denominator_ms,
      projector = denominator$H
    ),
    estimability = specification[c(
      "estimable", "fully_estimable", "df", "nominal_df",
      "aliased_columns", "reason", "projector_error"
    )],
    adjusted_for = specification$adjusted_for,
    containing_interactions = specification$containing_interactions,
    statistic = observed, p_value = p_value,
    exchangeable_units = "observations",
    permutation_strata = "unrestricted",
    permutation_algorithm = paste0(
      "term-specific Freedman-Lane permutation for Type ",
      specification$hypothesis_type, " hypothesis"
    ),
    attainable_permutations = plan$attainable,
    permutations_used = nrow(plan$orders),
    unique_permutations_used = if (nrow(plan$orders)) {
      nrow(unique(plan$orders))
    } else 0L,
    exact_permutation = plan$exact && nrow(plan$orders) > 0L,
    approximation_status = if (!nrow(plan$orders)) {
      "not tested"
    } else if (plan$exact) {
      "exhaustive term-specific residual permutation"
    } else {
      "Monte Carlo term-specific residual permutation"
    },
    warnings = notes,
    permutation_statistics = statistics
  )
  class(out) <- "permanova_test_spec"
  out
}

#' Unbalanced fixed-effect PERMANOVA
#'
#' Tests explicit estimable Type I, II, or III hypotheses over the final
#' full-model residual mean square. Type III requires a full-rank model and
#' explicit sum-to-zero contrasts when interactions are present. Type II
#' respects marginality: a term is adjusted for every term that does not
#' contain it and the higher-order interactions containing it are excluded, so
#' main effects remain estimable in the presence of an interaction (matching
#' the SAS/`car` "Type II" hypothesis).
#'
#' @param design Canonical fixed-effect Design Object.
#' @param hypothesis One of `"I"`, `"II"`, or `"III"`.
#' @param terms Optional fixed terms to test.
#' @param permutations A non-negative integer or `"all"`.
#' @param seed Optional permutation seed.
#' @param keep_permutations Retain permutation statistics.
#' @param max_exact Maximum exact permutation space.
#' @param tol Numerical tolerance.
#' @return An object of class `permanova_unbalanced`.
#' @export
permanova_unbalanced <- function(
    design, hypothesis = c("I", "II", "III"), terms = NULL,
    permutations = 999L, seed = NULL, keep_permutations = FALSE,
    max_exact = 100000L, tol = 1e-10) {
  hypothesis <- match.arg(hypothesis)
  specifications <- derive_unbalanced_hypotheses(
    design, hypothesis = hypothesis, terms = terms, tol = tol
  )
  .validate_unbalanced_hypotheses(specifications)
  full <- design$projections$fixed
  denominator <- list(
    H = diag(design$n) - full$H,
    df = design$n - full$rank
  )
  if (!is.null(seed)) set.seed(seed)
  plan <- .permutation_orders(design$n, permutations, max_exact)
  if (nrow(plan$orders) > 0L && is.finite(plan$attainable) &&
      plan$attainable < 100) {
    warning(
      "Only ", plan$attainable,
      " unrestricted observation permutations are attainable; interpret p-values cautiously.",
      call. = FALSE
    )
  }
  tests <- setNames(lapply(names(specifications), function(term) {
    .unbalanced_test_specification(
      term, specifications[[term]], design, denominator, plan, tol
    )
  }), names(specifications))
  if (!keep_permutations) {
    tests <- lapply(tests, function(test) {
      test$permutation_statistics <- NULL
      class(test) <- "permanova_test_spec"
      test
    })
  }
  table <- do.call(rbind, lapply(tests, function(test) {
    data.frame(
      term = test$term, hypothesis = test$hypothesis_type,
      df = test$numerator$df, nominal_df = test$estimability$nominal_df,
      ss = test$numerator$ss, ms = test$numerator$ms,
      denominator_df = test$denominator$df,
      pseudo_f = test$statistic, p_perm = test$p_value,
      stringsAsFactors = FALSE
    )
  }))
  table <- rbind(table, data.frame(
    term = "Residual", hypothesis = NA_character_,
    df = denominator$df, nominal_df = denominator$df,
    ss = trace_hgh(denominator$H, design$G),
    ms = trace_hgh(denominator$H, design$G) / denominator$df,
    denominator_df = NA_real_, pseudo_f = NA_real_, p_perm = NA_real_
  ))
  rownames(table) <- NULL
  out <- list(
    call = match.call(), design = design, hypothesis = hypothesis,
    table = table, tests = tests, hypotheses = specifications,
    estimability = list(
      fixed_rank = design$estimability$fixed_rank,
      fixed_columns = design$estimability$fixed_columns,
      aliased_columns = design$estimability$aliased,
      tested_terms = names(specifications)
    ),
    supported_design = paste0(
      "unbalanced fixed-effect model with fully estimable Type ",
      hypothesis, " hypotheses and a full-model residual denominator"
    ),
    permutation = list(
      algorithm = paste0(
        "term-specific Type ", hypothesis,
        " Freedman-Lane unrestricted residual permutation"
      ),
      exact = plan$exact && nrow(plan$orders) > 0L,
      attainable = plan$attainable,
      used = nrow(plan$orders), seed = seed
    )
  )
  class(out) <- "permanova_unbalanced"
  out
}

#' @export
print.permanova_unbalanced <- function(x, ...) {
  cat("Unbalanced fixed-effect PERMANOVA\n")
  cat("Hypothesis:", x$hypothesis, "\n")
  cat("Supported design:", x$supported_design, "\n")
  print(x$table, row.names = FALSE)
  invisible(x)
}
