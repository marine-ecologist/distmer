#' Continuous-covariate random-slope PERMANOVA
#'
#' Tests the population slope of a continuous covariate that also carries a
#' random slope across groups, i.e. `distance ~ x + (1 + x | group)` with `x`
#' numeric -- the covariate-by-random model of Anderson, Gorley & Clarke (2008,
#' Fig. 1.48b).
#'
#' The fixed slope has no exact single-source pseudo-F denominator in the full
#' observation space, because its expected mean square mixes the random-slope
#' and residual components with different coefficients. The test therefore
#' aggregates, exactly as the two-level categorical engine reduces to
#' per-subject difference vectors: it forms one within-group slope vector per
#' group (in the distance embedding), then tests whether the equally weighted
#' mean slope vector differs from zero against the between-group variation of
#' those slope vectors. Inference is by whole-group sign flips of the slope
#' vectors, valid under central symmetry of the group slopes under the null. A
#' parametric `F(1, ng - 1)` reference is added only for rank-one Euclidean
#' slope geometry with equal within-group covariate information (`SSx`), where
#' the statistic is the classical one-sample squared t statistic under normal,
#' identically distributed group slopes. Euclidean geometry alone is not enough
#' to justify that reference for multivariate or heteroscedastic slopes.
#'
#' The within-group slope is invariant to centring `x` within a group, so the
#' raw or a within-centred covariate ([decompose_covariate()]) give the same
#' test; the covariate need only vary inside each group.
#'
#' @param design A compiled Design Object with exactly one `(1 + x | group)`
#'   random term whose slope variable is numeric.
#' @param permutations A non-negative integer or `"all"` (whole-group sign
#'   flips).
#' @param seed Optional sign-flip seed.
#' @param keep_permutations Retain the sign-flip statistics.
#' @param max_exact Maximum exact sign-flip space.
#' @param tol Numerical tolerance.
#' @return An object of class `permanova_random_slope_continuous`.
#' @export
permanova_random_slope_continuous <- function(
    design, permutations = 999L, seed = NULL, keep_permutations = FALSE,
    max_exact = 100000L, tol = 1e-10) {
  if (!inherits(design, "permanova_design") || is.null(design$G)) {
    stop("`design` must be a compiled Design Object containing `distance`.",
         call. = FALSE)
  }
  if (length(design$random_terms) != 1L ||
      design$random_terms[[1L]]$origin != "slope") {
    stop("Continuous random-slope inference requires one `(1 + x | group)` term.",
         call. = FALSE)
  }
  spec <- design$random_terms[[1L]]
  slope_name <- spec$slope
  group_name <- spec$group
  x <- design$data[[slope_name]]
  if (!is.numeric(x) || is.factor(x)) {
    stop("The random-slope variable `", slope_name,
         "` must be numeric; a two-level factor uses `permanova_random_slope()`.",
         call. = FALSE)
  }
  fixed_variables <- all.vars(design$fixed_formula)
  if (length(design$fixed_effects) != 1L || length(fixed_variables) != 1L ||
      !identical(fixed_variables[[1L]], slope_name)) {
    stop("The continuous random slope requires its covariate as the sole fixed effect.",
         call. = FALSE)
  }
  group <- droplevels(factor(design$data[[group_name]]))
  groups <- levels(group)
  ng <- length(groups)
  if (ng < 3L) {
    stop("At least three groups are required to estimate a random slope.",
         call. = FALSE)
  }
  x <- as.numeric(x)
  n <- design$n
  G <- design$G

  # One within-group slope contrast per group: w_g[i] = (x_i - xbar_g) / SS_x(g)
  # for observations i in group g. Because the weights sum to zero within the
  # group, applying them to the embedding removes the group mean automatically,
  # so W %*% Y is the vector of within-group OLS slopes.
  W <- matrix(0, nrow = ng, ncol = n)
  ssx_by_group <- numeric(ng)
  for (k in seq_len(ng)) {
    idx <- which(group == groups[[k]])
    xc <- x[idx] - mean(x[idx])
    ssx <- sum(xc^2)
    if (ssx <= tol) {
      stop("Group `", groups[[k]], "` has no within-group covariate variation; ",
           "its slope is not estimable.", call. = FALSE)
    }
    ssx_by_group[[k]] <- ssx
    W[k, idx] <- xc / ssx
  }
  K <- W %*% G %*% t(W)
  K <- (K + t(K)) / 2

  H_mean <- matrix(1 / ng, ng, ng)
  H_between <- diag(ng) - H_mean
  numerator_ss <- trace_hgh(H_mean, K)
  denominator_ss <- trace_hgh(H_between, K)
  denominator_df <- ng - 1L
  denominator_ms <- denominator_ss / denominator_df
  if (!is.finite(denominator_ms) || abs(denominator_ms) <= tol) {
    stop("The between-group random-slope denominator is numerically zero.",
         call. = FALSE)
  }
  observed <- numerator_ss / denominator_ms

  if (!is.null(seed)) set.seed(seed)
  plan <- .sign_flip_plan(ng, permutations, max_exact)
  statistics <- vapply(seq_len(nrow(plan$signs)), function(index) {
    signs <- plan$signs[index, ]
    permuted <- K * tcrossprod(signs)
    den <- trace_hgh(H_between, permuted) / denominator_df
    if (!is.finite(den) || abs(den) <= tol) {
      stop("A sign-flipped random-slope denominator is numerically zero.",
           call. = FALSE)
    }
    trace_hgh(H_mean, permuted) / den
  }, numeric(1))
  p_value <- NA_real_
  if (length(statistics)) {
    extreme <- sum(statistics >= observed - tol, na.rm = TRUE)
    p_value <- if (plan$exact) extreme / length(statistics) else
      (extreme + 1) / (length(statistics) + 1)
  }
  slope_eigen <- eigen(K, symmetric = TRUE, only.values = TRUE)$values
  slope_scale <- max(abs(slope_eigen))
  slope_rank <- if (is.finite(slope_scale) && slope_scale > 0) {
    sum(slope_eigen > max(100 * tol, sqrt(.Machine$double.eps)) * slope_scale)
  } else 0L
  information_error <- max(abs(ssx_by_group - mean(ssx_by_group))) /
    mean(ssx_by_group)
  equal_information <- information_error <=
    max(100 * tol, sqrt(.Machine$double.eps))
  parametric_supported <- isTRUE(design$spectral$euclidean) &&
    slope_rank == 1L && equal_information
  parametric_p <- if (parametric_supported) {
    stats::pf(observed, 1L, denominator_df, lower.tail = FALSE)
  } else NA_real_

  # Between-group slope mean square minus the sampling contribution gives a
  # method-of-moments slope variance component; the residual mean square is the
  # pooled within-group residual after both the group means and slopes.
  Zg <- stats::model.matrix(~ 0 + group)
  H_full <- projection_matrix(cbind(1, x, Zg, Zg * x), tol = tol)$H
  residual_df <- n - as.integer(round(sum(diag(H_full))))
  residual_ms <- if (residual_df > 0L) {
    trace_hgh(diag(n) - H_full, G) / residual_df
  } else NA_real_
  slope_component <- denominator_ms - if (is.finite(residual_ms)) {
    residual_ms * mean(1 / ssx_by_group)
  } else 0

  notes <- paste(
    "Whole-group sign flips of within-group slope vectors require central",
    "symmetry of the group slopes under the null population slope of zero."
  )
  if (is.finite(plan$attainable) && plan$attainable < 100) {
    notes <- c(notes, paste0("Only ", plan$attainable,
                             " group sign assignments are attainable."))
  }
  if (!isTRUE(design$spectral$euclidean)) {
    notes <- c(notes, paste(
      "The complete non-Euclidean Gower spectrum was retained; the parametric",
      "p-value was suppressed."))
  }
  if (isTRUE(design$spectral$euclidean) && slope_rank != 1L) {
    notes <- c(notes, paste(
      "The Euclidean slope geometry is not rank one; the univariate F(1, ng-1)",
      "reference was suppressed."))
  }
  if (isTRUE(design$spectral$euclidean) && !equal_information) {
    notes <- c(notes, paste(
      "Within-group covariate SSx differs among groups; the equal-information",
      "F(1, ng-1) reference was suppressed."))
  }

  test <- list(
    term = slope_name, kind = "continuous random slope",
    hypothesis = paste0("the population slope of `", slope_name, "` is zero"),
    estimand = "equally weighted within-group covariate slopes",
    numerator = list(ss = numerator_ss, df = 1L, ms = numerator_ss,
                     projector = H_mean),
    denominator = list(
      expression = "MS[between-group random slopes]",
      ss = denominator_ss, df = denominator_df, ms = denominator_ms,
      projector = H_between),
    statistic = observed, p_value = p_value, p_value_parametric = parametric_p,
    exchangeable_units = group_name,
    permutation_algorithm = "whole-group sign flips of within-group slope vectors",
    exchangeability_assumption = "central symmetry of group slopes under the null",
    parametric_assumption = if (parametric_supported) {
      "normal identically distributed rank-one group slopes with equal within-group SSx"
    } else "not satisfied; parametric p-value suppressed",
    attainable_permutations = plan$attainable,
    permutations_used = nrow(plan$signs),
    exact_permutation = plan$exact && nrow(plan$signs) > 0L,
    warnings = notes,
    permutation_statistics = if (keep_permutations) statistics else NULL
  )
  class(test) <- "permanova_test_spec"

  table <- data.frame(
    term = c(slope_name, "Between-group random slope", "Residual"),
    df = c(1L, denominator_df, residual_df),
    ss = c(numerator_ss, denominator_ss,
           if (is.finite(residual_ms)) residual_ms * residual_df else NA_real_),
    ms = c(numerator_ss, denominator_ms, residual_ms),
    denominator_df = c(denominator_df, NA_real_, NA_real_),
    pseudo_f = c(observed, NA_real_, NA_real_),
    p_perm = c(p_value, NA_real_, NA_real_),
    p_parametric = c(parametric_p, NA_real_, NA_real_),
    stringsAsFactors = FALSE
  )

  out <- list(
    call = match.call(), design = design, table = table, test = test,
    tests = stats::setNames(list(test), slope_name),
    variance_components = list(
      unrestricted = c(slope = slope_component, residual = residual_ms),
      method = paste(
        "between-group slope mean square less residual MS times mean reciprocal",
        "within-group covariate SSx"
      ),
      within_group_ssx = stats::setNames(ssx_by_group, groups)
    ),
    estimand = test$estimand,
    supported_design = paste(
      "continuous covariate with a correlated random slope across groups;",
      "aggregated within-group slope vectors tested by whole-group sign flips"
    ),
    permutation = list(
      algorithm = test$permutation_algorithm,
      assumption = test$exchangeability_assumption,
      exact = test$exact_permutation, attainable = plan$attainable,
      used = nrow(plan$signs), seed = seed)
  )
  class(out) <- "permanova_random_slope_continuous"
  out
}

#' @export
print.permanova_random_slope_continuous <- function(x, ...) {
  cat("Continuous-covariate random-slope PERMANOVA\n")
  cat("Estimand:", x$estimand, "\n")
  cat("Denominator:", x$test$denominator$expression,
      sprintf("(df = %d)\n", x$test$denominator$df))
  print(x$table, row.names = FALSE)
  invisible(x)
}
