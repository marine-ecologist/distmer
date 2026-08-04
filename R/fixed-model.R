.prepare_fixed_model <- function(formula, data, tol = 1e-10) {
  if (!inherits(formula, "formula")) {
    stop("`formula` must be a formula.", call. = FALSE)
  }

  terms_object <- terms(formula, data = data)
  if (attr(terms_object, "response") == 1L) {
    terms_object <- delete.response(terms_object)
  }

  model_frame <- model.frame(
    terms_object,
    data = data,
    na.action = na.fail,
    drop.unused.levels = TRUE
  )
  X <- model.matrix(terms_object, data = model_frame)
  assign <- attr(X, "assign")
  term_labels <- attr(terms_object, "term.labels")

  intercept_cols <- which(assign == 0L)
  base_projection <- projection_matrix(X[, intercept_cols, drop = FALSE], tol = tol)

  cumulative <- vector("list", length(term_labels))
  previous <- base_projection

  for (k in seq_along(term_labels)) {
    cols <- which(assign %in% c(0L, seq_len(k)))
    current <- projection_matrix(X[, cols, drop = FALSE], tol = tol)
    cumulative[[k]] <- list(
      term = term_labels[[k]],
      H_term = current$H - previous$H,
      H_full = current$H,
      df_term = current$rank - previous$rank,
      rank_full = current$rank
    )
    previous <- current
  }

  list(
    terms = terms_object,
    model_frame = model_frame,
    X = X,
    term_labels = term_labels,
    base = base_projection,
    cumulative = cumulative,
    full = previous
  )
}

.fixed_sequential_table <- function(G, prepared) {
  n <- nrow(G)
  term_rows <- lapply(prepared$cumulative, function(term) {
    data.frame(
      term = term$term,
      df = term$df_term,
      ss = trace_hgh(term$H_term, G),
      stringsAsFactors = FALSE
    )
  })

  table <- do.call(rbind, term_rows)
  H_residual <- diag(n) - prepared$full$H
  residual_ss <- trace_hgh(H_residual, G)
  residual_df <- n - prepared$full$rank
  residual_ms <- residual_ss / residual_df

  table$ms <- table$ss / table$df
  table$pseudo_f <- table$ms / residual_ms
  table$p_perm <- NA_real_

  rbind(
    table,
    data.frame(
      term = "Residual",
      df = residual_df,
      ss = residual_ss,
      ms = residual_ms,
      pseudo_f = NA_real_,
      p_perm = NA_real_
    ),
    data.frame(
      term = "Total",
      df = n - 1L,
      ss = sum(diag(G)),
      ms = NA_real_,
      pseudo_f = NA_real_,
      p_perm = NA_real_
    )
  )
}

#' Fixed-effect sequential PERMANOVA
#'
#' Initial clean-room multi-factor engine. It implements Type-I/sequential
#' sums of squares and tests each term over the final residual mean square.
#' Permutations are unrestricted raw-data permutations represented by the
#' same row-and-column permutation of the Gower-centred matrix.
#'
#' This is not yet a replacement for the legacy program's EMS denominator
#' logic for mixed, random, or nested designs.
#'
#' @param distance A `dist` object or full distance matrix.
#' @param formula One-sided or two-sided model formula; the response is ignored.
#' @param data Model data in the same row order as the distance matrix.
#' @param permutations Number of random raw-data permutations.
#' @param seed Optional random seed.
#' @param keep_permutations Retain the matrix of permutation pseudo-F values.
#' @param tol QR tolerance.
#' @return An object of class `permanova_fixed`.
#' @export
permanova_fixed <- function(
    distance,
    formula,
    data,
    permutations = 999L,
    seed = NULL,
    keep_permutations = FALSE,
    tol = 1e-10) {
  distance <- as_distance_matrix(distance)
  n <- nrow(distance)

  if (nrow(data) != n) {
    stop("`data` rows must match the distance-matrix observation order.", call. = FALSE)
  }

  if (length(permutations) != 1L || permutations < 0 || permutations != as.integer(permutations)) {
    stop("`permutations` must be one non-negative integer.", call. = FALSE)
  }

  prepared <- .prepare_fixed_model(formula, data, tol = tol)
  if (nrow(prepared$model_frame) != n) {
    stop("The model frame does not contain exactly one row per distance-matrix observation.", call. = FALSE)
  }

  G <- gower_center(distance)
  observed_table <- .fixed_sequential_table(G, prepared)
  term_count <- length(prepared$term_labels)
  permutation_f <- matrix(
    NA_real_,
    nrow = as.integer(permutations),
    ncol = term_count,
    dimnames = list(NULL, prepared$term_labels)
  )

  if (permutations > 0L) {
    if (!is.null(seed)) {
      set.seed(seed)
    }

    progress <- .permutation_progress(
      as.integer(permutations), "Permuting observations"
    )
    on.exit(progress$close(), add = TRUE)
    for (b in seq_len(permutations)) {
      idx <- sample.int(n, size = n, replace = FALSE)
      G_perm <- G[idx, idx, drop = FALSE]
      perm_table <- .fixed_sequential_table(G_perm, prepared)
      permutation_f[b, ] <- perm_table$pseudo_f[seq_len(term_count)]
      progress$tick(b)
    }

    observed_f <- observed_table$pseudo_f[seq_len(term_count)]
    observed_table$p_perm[seq_len(term_count)] <- vapply(
      seq_len(term_count),
      function(k) {
        (1 + sum(permutation_f[, k] >= observed_f[[k]])) /
          (as.integer(permutations) + 1)
      },
      numeric(1)
    )
  }

  out <- list(
    call = match.call(),
    table = observed_table,
    formula = prepared$terms,
    permutations = as.integer(permutations),
    seed = seed,
    method = "unrestricted raw-data permutation; sequential SS; residual denominator",
    design_rank = prepared$full$rank,
    permutation_f = if (keep_permutations) permutation_f else NULL
  )
  class(out) <- "permanova_fixed"
  out
}

#' @export
print.permanova_fixed <- function(x, ...) {
  cat("Fixed-effect sequential PERMANOVA\n")
  cat("Method:", x$method, "\n")
  print(x$table, row.names = FALSE)
  invisible(x)
}

.fix2_distinct_perms <- function(Fperm, F_obs = NULL, tol = 1e-9) {
  vals <- c(F_obs, Fperm)
  vals <- vals[is.finite(vals)]
  if (!length(vals)) return(NA_integer_)
  # bucket to tolerance so floating-point equal statistics count as one
  length(unique(round(vals / max(abs(vals)), digits = ceiling(-log10(tol)))))
}

# estimable df within-stratum is 0.
.fix3_confounded <- function(term_values, stratum_values) {
  # for each stratum unit, does the term vary?
  by_unit <- tapply(as.character(term_values), stratum_values,
                    function(v) length(unique(v)))
  all(by_unit == 1L)                        # constant in every unit -> confounded
}

# When confounded, emit a 0-df / no-test row instead of fitting:
.fix3_no_test_row <- function(term_name) {
  data.frame(
    Source   = term_name,
    df       = 0L,
    SS       = 0,
    MS       = NA_real_,
    `Pseudo-F` = NA_real_,
    `P(perm)`  = NA_real_,
    note     = "No test (confounded with random stratum)",
    check.names = FALSE, stringsAsFactors = FALSE
  )
}


# Classify a random grouping relative to a fixed factor:
#   "whole_plot"  : fixed factor constant within each group (between-unit)
#   "block"       : fixed factor varies within each group (within-unit / RCB)
#   "partial"     : mixed -> genuine confounding (keep the refusal, that's D2/B2)
.fix4_classify_grouping <- function(fixed_values, group_values) {
  by_grp <- tapply(as.character(fixed_values), group_values,
                   function(v) length(unique(v)))
  n_lev  <- length(unique(as.character(fixed_values)))
  if (all(by_grp == 1L))            return("whole_plot")   # constant in every group
  if (all(by_grp == n_lev))         return("block")        # full factor in every group
  "partial"                                                # some groups partial -> confounded
}


.fix5_no_fixed_message <- function(random_terms) {
  .stopf(paste0(
    "No fixed factor in the model. distmer tests fixed effects against random ",
    "strata; it does not estimate a pure random-effects (variance-components) ",
    "decomposition.\n",
    "  Your random structure: %s\n",
    "  To analyse this design, either:\n",
    "   - add the fixed factor of interest, e.g. `d ~ treatment + (1 | %s)`, or\n",
    "   - for pure variance components, use a dedicated nested-ANOVA / REML tool.\n",
    "  (See the 'designs distmer does not fit' vignette section.)"),
    paste(random_terms, collapse = ", "),
    random_terms[[1]] %||% "group")
}
