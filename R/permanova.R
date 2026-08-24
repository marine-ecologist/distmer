.or <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

# TRUE for the simplest single-random-intercept design: exactly one fixed factor
# (categorical, not a covariate) that is constant within every grouping unit --
# the classic whole-plot test handled by the exact-EMS `permanova_mixed()`.
# Anything richer (more fixed terms, a within-unit factor, or a covariate) is a
# split-plot / RCB / ANCOVA and routes to `permanova_mixed_covariate()`.
.is_simple_whole_plot <- function(design) {
  labels <- attr(stats::terms(design$fixed_formula), "term.labels")
  variables <- all.vars(design$fixed_formula)
  if (length(labels) != 1L || length(variables) != 1L) return(FALSE)
  values <- design$data[[variables[[1L]]]]
  if (is.numeric(values) && !is.factor(values)) return(FALSE)
  unit <- design$groups[[design$random_terms[[1L]]$label]]
  all(tapply(as.character(values), unit,
             function(v) length(unique(v))) == 1L)
}

# Route a compiled Design Object to the engine that matches its provenance.
# Returns the engine label, the engine function, and any engine-specific fixed
# arguments. Designs outside every supported envelope stop with guidance.
.permanova_dispatch <- function(design, type) {
  origins <- if (length(design$random_terms)) {
    vapply(design$random_terms, `[[`, character(1), "origin")
  } else character(0)

  # A variable used as both a fixed effect and a random *grouping* factor is a
  # misspecification -- almost always `(1 | fixedfactor / unit)`, where the `/`
  # was meant to nest `unit`, not to make the fixed factor random. (The random
  # slope's own covariate is fine: it is not a grouping variable.)
  if (length(origins)) {
    # Only `/`-style main/nesting groupings are flagged: a fixed factor that is a
    # component of a random *interaction* (e.g. `(1 | unit:fixedfactor)`) is a
    # legitimate random interaction, not a misdeclared grouping.
    grouping_vars <- unique(unlist(lapply(design$random_terms, function(t) {
      if (isTRUE(t$origin %in% c("main", "nesting"))) t$variables else NULL
    })))
    conflict <- intersect(all.vars(design$fixed_formula), grouping_vars)
    if (length(conflict)) {
      stop(
        "The factor(s) ", paste(conflict, collapse = ", "), " are declared as ",
        "both fixed effects and random grouping factors. In an `lmer`-style ",
        "formula `(1 | a / b)` makes BOTH `a` and `b` random. If these are fixed ",
        "factors and only the innermost unit is a random effect nested within ",
        "their combination, drop the `/` and use one `(1 | unit)` term with unit ",
        "labels unique to each factor combination, e.g. `d ~ a * b + (1 | unit)`.",
        call. = FALSE
      )
    }
  }

  # distmer tests fixed effects against random strata. With no fixed term there
  # is nothing to test, and the engines would otherwise return an SS table with
  # no statistic -- which reads like a result. Refuse, and name the function
  # that does decompose a pure random structure.
  if (length(origins) && !length(design$fixed_effects)) {
    nested <- length(origins) > 1L && all(origins %in% c("main", "nesting")) &&
      sum(origins == "main") == 1L
    estimator <- if (nested) {
      "estimate_nested_variance_components()"
    } else "estimate_variance_components()"
    stop(
      "No fixed term in the model. distmer tests fixed effects against random ",
      "strata; it does not report a pure random-effects decomposition through ",
      "`permanova()`. Either add the fixed factor of interest, for example ",
      "`d ~ treatment + ",
      paste0("(1 | ", design$random_terms[[1L]]$label, ")"),
      "`, or call `", estimator, "` on the compiled design for the ",
      "variance components themselves.",
      call. = FALSE
    )
  }

  if (!length(origins)) {
    return(list(engine = "fixed effects", fn = permanova_unbalanced,
                fixed_args = list(hypothesis = .or(type, "I"))))
  }
  if (length(origins) == 1L && origins == "slope") {
    slope_values <- design$data[[design$random_terms[[1L]]$slope]]
    if (is.numeric(slope_values) && !is.factor(slope_values)) {
      return(list(engine = "continuous random slope",
                  fn = permanova_random_slope_continuous, fixed_args = list()))
    }
    return(list(engine = "random slope", fn = permanova_random_slope,
                fixed_args = list()))
  }
  if (length(origins) == 1L && origins == "main") {
    if (all(design$balance$random_balanced)) {
      # The exact single-whole-plot random-intercept engine keeps its richer
      # symbolic-EMS output; any other balanced single-random-intercept design
      # (extra fixed factors, a within-unit/sub-plot factor, or a covariate) goes
      # to the general split-plot / RCB / ANCOVA engine.
      if (.is_simple_whole_plot(design)) {
        return(list(engine = "balanced random intercept", fn = permanova_mixed,
                    fixed_args = list()))
      }
      return(list(engine = "balanced split-plot mixed",
                  fn = permanova_mixed_covariate,
                  fixed_args = list(hypothesis = .or(type, "I"))))
    }
    # Unbalanced single random intercept: a lone whole-plot factor keeps the
    # unit-centroid engine; anything richer (extra fixed factors, a within-unit
    # sub-plot factor, or a covariate) goes to the unbalanced split-plot engine.
    if (.is_simple_whole_plot(design)) {
      return(list(engine = "unequal random intercept",
                  fn = permanova_unequal_random_intercept, fixed_args = list()))
    }
    return(list(engine = "unbalanced split-plot mixed",
                fn = permanova_unbalanced_splitplot,
                fixed_args = list(hypothesis = .or(type, "I"))))
  }
  if (all(origins %in% c("main", "nesting")) && sum(origins == "main") == 1L) {
    has_covariate <- any(vapply(all.vars(design$fixed_formula), function(variable) {
      value <- design$data[[variable]]
      is.numeric(value) && !is.factor(value)
    }, logical(1)))
    if (all(design$balance$random_balanced)) {
      if (has_covariate) {
        return(list(engine = "balanced nested with covariate",
                    fn = permanova_nested_covariate,
                    fixed_args = list(hypothesis = .or(type, "I"))))
      }
      return(list(engine = "balanced nested", fn = permanova_nested,
                  fixed_args = list()))
    }
    # Deeper than parent/child: cascaded equal-weight centroids at the
    # treatment's own hierarchy level. The two-level engines keep their
    # validated observation-weighted parent centroids.
    if (length(design$random_terms) > 2L) {
      return(list(engine = "unbalanced nested (multi-level)",
                  fn = permanova_unbalanced_nested_multi,
                  fixed_args = list(hypothesis = .or(type, "I"))))
    }
    assignment <- .unbalanced_nested_assignment(design)
    if (assignment$level == "top") {
      return(list(engine = "unbalanced nested (top unit)",
                  fn = permanova_unbalanced_nested, fixed_args = list()))
    }
    return(list(
      engine = paste0("unbalanced nested (", assignment$level, " level)"),
      fn = permanova_unbalanced_nested_lower,
      fixed_args = list(level = assignment$level)
    ))
  }
  if (any(origins == "nesting") && any(origins == "interaction")) {
    if (all(design$balance$random_balanced)) {
      return(list(engine = "balanced crossed-nested hybrid",
                  fn = permanova_hybrid, fixed_args = list()))
    }
    return(list(engine = "unbalanced crossed-nested hybrid",
                fn = permanova_unbalanced_hybrid, fixed_args = list()))
  }
  # `(1 | unit) + (1 | unit:subplot)`: a unit random effect plus its interaction
  # with a sub-plot factor -- the asymmetric nested mixed design (three strata).
  if (length(origins) == 2L && sum(origins == "main") == 1L &&
      sum(origins == "interaction") == 1L) {
    main_var <- design$random_terms[[which(origins == "main")]]$variables
    inter_vars <- design$random_terms[[which(origins == "interaction")]]$variables
    if (length(main_var) == 1L && length(inter_vars) > 1L &&
        main_var %in% inter_vars) {
      return(list(engine = "asymmetric nested mixed",
                  fn = permanova_asymmetric_mixed,
                  fixed_args = list(hypothesis = .or(type, "I"))))
    }
  }
  # A `/`-generated chain reaching this point carries crossed random terms that
  # do not form a validated region/site hybrid. The crossed engine would take
  # the route and then refuse from inside, so refuse here where the structure is
  # still visible and can be named.
  if (any(origins == "nesting")) {
    stop(
      "This design mixes a `/`-generated nesting chain (",
      paste(design$random_terms[[which(origins == "nesting")[[1L]]]]$label,
            collapse = ", "),
      ") with crossed random terms. distmer does not infer a denominator for a ",
      "nesting chain crossed with an unrelated random factor. Either drop the ",
      "crossed term and use the nested engine, or declare the validated ",
      "crossed-nested hybrid `(1 | region/site) + (1 | treatment:region)`.",
      call. = FALSE
    )
  }
  if (any(origins == "interaction") || sum(origins == "main") >= 2L) {
    # `permanova_crossed()` covers the complete, equally replicated case and
    # refuses anything else. When a cell is missing or replication is uneven,
    # route to the engine built for that instead of surfacing its refusal: the
    # two agree exactly wherever both apply.
    if (.crossed_is_regular(design)) {
      return(list(engine = "crossed / randomized block", fn = permanova_crossed,
                  fixed_args = list(p_method = "permutation")))
    }
    return(list(engine = "crossed mixed (incomplete / unequally replicated)",
                fn = permanova_crossed_incomplete, fixed_args = list()))
  }
  stop(
    "permanova(): this random-effect structure is not auto-routed yet. ",
    "Compile the design and call a specific engine (see `?permanova`).",
    call. = FALSE
  )
}

# Is every cell of the crossed factorial present, with the same number of
# observations in each? That is what `permanova_crossed()` requires, so it
# decides which crossed engine a design is routed to.
.crossed_is_regular <- function(design) {
  fixed_variables <- all.vars(design$fixed_formula)
  random_variables <- unique(unlist(lapply(design$random_terms, `[[`, "variables")))
  names <- unique(c(fixed_variables, random_variables))
  if (!length(names) || !all(names %in% names(design$data))) return(FALSE)
  columns <- design$data[names]
  if (!all(vapply(columns, is.factor, logical(1)))) return(FALSE)
  counts <- table(interaction(lapply(columns, droplevels), drop = TRUE))
  expected <- prod(vapply(columns, function(z) nlevels(droplevels(z)), integer(1)))
  length(counts) == expected && length(unique(as.integer(counts))) == 1L
}

# Set sum-to-zero contrasts on the factors named in a formula, on a copy of the
# data. This makes Type III main-effect hypotheses (and cross-level factor
# interactions) the correct marginal tests without the caller touching
# `contrasts()` or session-wide options. Type I/II factor sums of squares are
# contrast-invariant, so this never changes them; grouping factors use their
# own indicator coding downstream and are unaffected. Numeric covariates are
# left untouched.
.apply_sum_contrasts <- function(formula, data) {
  if (is.null(data)) return(data)
  data <- as.data.frame(data)
  variables <- intersect(all.vars(formula), names(data))
  for (variable in variables) {
    column <- data[[variable]]
    if (is.factor(column)) {
      column <- droplevels(column)
      if (nlevels(column) >= 2L) {
        stats::contrasts(column) <- stats::contr.sum(nlevels(column))
        data[[variable]] <- column
      }
    }
  }
  data
}

#' Fit a distance-based PERMANOVA for any supported design
#'
#' A single front door that compiles a formula into a canonical Design Object
#' (or accepts one directly), routes it to the engine matching its
#' random-effect provenance and balance, and returns a result that prints as a
#' PERMANOVA table. No denominator or permutation scheme is ever chosen
#' silently: a design outside every supported envelope stops with a message
#' naming the reason.
#'
#' @param formula A model formula, or a compiled `permanova_design` object. The
#'   dissimilarity may be given on the left-hand side, as in
#'   `bray ~ treatment + (1 | site)` (the idiom used by `vegan::adonis2()`), or
#'   the formula may be one-sided (`~ treatment + (1 | site)`) with the
#'   dissimilarity passed through `distance`.
#' @param data A data frame with one row per observation. Ignored when
#'   `formula` is already a Design Object.
#' @param distance A `dist` object or distance matrix. Needed only when the
#'   dissimilarity is not given on the formula's left-hand side.
#' @param type Sum-of-squares partition for the fixed-effect terms: `"I"`
#'   (sequential), `"II"` (conditional), or `"III"` (partial). Applies to
#'   fixed-effect designs; for designs whose denominator is fixed by the random
#'   structure it is reported as not applicable.
#' @param permutations A non-negative integer or `"all"`.
#' @param mc_permutations Number of Monte-Carlo draws for the asymptotic
#'   `P(MC)` p-value, reported alongside `P(perm)` for designs tested in the
#'   full observation space. `P(MC)` stays smooth when the true permutation
#'   space is small (the usual case for the upper strata of a mixed design).
#'   Set to `0` to disable. Reproducible from `seed`.
#' @param seed Optional random seed. Also seeds the `P(MC)` draws.
#' @param keep_permutations Retain permutation statistics on the fitted object.
#' @param truncate Also report zero-truncated variance components where the
#'   engine estimates them.
#' @param max_exact Maximum permutation space allowed for exact enumeration.
#' @param tol Numerical tolerance.
#' @return An object of class `permanova` wrapping the engine result. Use
#'   `summary()` for the full table with variance components.
#' @examples
#' data(dune, dune.env, package = "vegan")
#' bray <- vegan::vegdist(dune, method = "bray")
#' permanova(bray ~ Management, dune.env, permutations = 99, seed = 1)
#' @export
permanova <- function(formula, data = NULL, distance = NULL, type = NULL,
                      permutations = 999L, mc_permutations = 9999L, seed = NULL,
                      keep_permutations = FALSE, truncate = FALSE,
                      max_exact = 100000L, tol = 1e-10) {
  if (inherits(formula, "permanova_design")) {
    design <- formula
  } else if (inherits(formula, "formula")) {
    model_formula <- formula
    # `distance ~ terms`: the left-hand side supplies the dissimilarity, as in
    # vegan::adonis2(). Evaluate it (in `data`, then the formula's environment)
    # and strip it, leaving a one-sided model formula. A separately supplied
    # `distance` still works with the one-sided `~ terms` form.
    if (length(formula) == 3L && is.null(distance)) {
      distance <- eval(formula[[2L]], envir = data,
                       enclos = environment(formula))
      model_formula <- stats::as.formula(
        call("~", formula[[3L]]), env = environment(formula)
      )
    }
    if (is.null(distance)) {
      stop(
        "Supply the dissimilarity: `permanova(D ~ terms, data)` or ",
        "`permanova(~ terms, data, distance = D)`.", call. = FALSE
      )
    }
    data <- .apply_sum_contrasts(model_formula, data)
    design <- compile_permanova_design(model_formula, data, distance, tol = tol)
  } else {
    stop("`formula` must be a formula or a `permanova_design` object.",
         call. = FALSE)
  }
  if (!is.null(type)) type <- match.arg(type, c("I", "II", "III"))

  route <- .permanova_dispatch(design, type)
  options_list <- c(
    list(permutations = permutations, seed = seed,
         keep_permutations = keep_permutations, truncate = truncate,
         max_exact = max_exact, tol = tol),
    route$fixed_args
  )
  usable <- intersect(names(options_list), names(formals(route$fn)))
  fit <- do.call(route$fn, c(list(design), options_list[usable]))

  applies_type <- "hypothesis" %in% usable
  if (!is.null(type) && !applies_type) {
    warning(
      "`type` controls fixed-effect partitioning; it does not apply to the ",
      route$engine, " design, whose denominator is set by the random ",
      "structure. It was ignored.", call. = FALSE
    )
  }
  out <- list(
    call = match.call(), design = design, engine = route$engine,
    type = if (applies_type) .or(type, "I") else NA_character_,
    fit = fit, seed = seed, permutations = permutations,
    mc_permutations = mc_permutations, keep_permutations = keep_permutations,
    max_exact = max_exact, tol = tol
  )
  class(out) <- "permanova"
  # Convenience: `$table` is the assembled PERMANOVA table (Source/df/SS/MS/
  # Pseudo-F/P(perm)/Unique perms); the raw engine table stays at `$fit$table`.
  out$table <- .permanova_table(out)
  out
}

.permanova_tests <- function(fit) {
  if (!is.null(fit$tests)) return(fit$tests)
  if (!is.null(fit$test)) return(stats::setNames(list(fit$test), fit$test$term))
  list()
}

.permanova_used <- function(fit) {
  .or(fit$permutation$used,
      .or(fit$permutation$attainable,
          .or(fit$test$permutations_used, NA_real_)))
}

# Monte-Carlo asymptotic p-value P(MC) for a distance-based pseudo-F.
# Under exchangeability the numerator and denominator sums of squares are each
# asymptotically distributed as a linear combination of independent chi-squares
# weighted by the eigenvalues of the Gower-centred matrix G (McArdle & Anderson
# 2001; Anderson & Robinson 2003). Sampling that theoretical distribution gives
# a p-value that stays smooth when the true permutation space is tiny -- the
# usual case for the upper strata of a mixed design. With a single non-zero
# eigenvalue (a univariate Euclidean response) it reduces exactly to the
# classical F distribution: numerator and denominator become the same lambda
# times independent chi-squares on df_num and df_den, so the ratio is F.
.monte_carlo_p <- function(lambda, df_num, df_den, f_obs, nmc = 9999L,
                           tol = 1e-9) {
  if (!length(lambda) || !is.finite(f_obs) || !is.finite(df_num) ||
      !is.finite(df_den) || df_num < 1 || df_den < 1) return(NA_real_)
  r <- length(lambda)
  num <- as.numeric(lambda %*% matrix(stats::rchisq(r * nmc, df_num), r, nmc)) / df_num
  den <- as.numeric(lambda %*% matrix(stats::rchisq(r * nmc, df_den), r, nmc)) / df_den
  fmc <- num / den
  keep <- is.finite(fmc)
  if (!any(keep)) return(NA_real_)
  (sum(fmc[keep] >= f_obs - max(100 * tol, sqrt(.Machine$double.eps)) * max(1, abs(f_obs))) + 1) / (sum(keep) + 1)
}

# P(MC) for every row of the assembled table, seeded reproducibly from the
# object. Only the full-observation-space designs (those carrying a Total row)
# are eligible; aggregated-unit engines test in a reduced centroid space where
# the full-G eigenvalue asymptotics do not apply, so their P(MC) is left NA.
.permanova_mc_column <- function(object, source, df_num, df_den, f_obs, enabled) {
  n_mc <- .or(object$mc_permutations, 0L)
  out <- rep(NA_real_, length(f_obs))
  if (!isTRUE(enabled) || !is.finite(n_mc) || n_mc <= 0L) return(out)
  tests <- .permanova_tests(object$fit)
  old <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv)
  } else NULL
  on.exit(if (is.null(old)) {
    suppressWarnings(rm(".Random.seed", envir = .GlobalEnv))
  } else assign(".Random.seed", old, envir = .GlobalEnv))
  set.seed(.or(object$seed, 1L))
  for (i in seq_along(f_obs)) {
    if (is.finite(f_obs[i]) && is.finite(df_den[i]) && is.finite(df_num[i])) {
      spec <- tests[[source[[i]]]]
      if (is.null(spec) && length(tests)) {
        match_term <- which(vapply(tests, function(x) {
          identical(as.character(x$term), source[[i]])
        }, logical(1)))
        if (length(match_term) == 1L) spec <- tests[[match_term]]
      }
      H_num <- spec$numerator$projector
      H_den <- spec$denominator$projector
      n <- object$design$n
      if (is.null(H_num) || is.null(H_den) ||
          !all(dim(H_num) == c(n, n)) || !all(dim(H_den) == c(n, n))) next
      # The asymptotic spectrum must be restricted to the exact numerator plus
      # denominator stratum. Using the full observation-level G here lets
      # irrelevant lower-stratum variation change an upper-stratum p-value.
      H_test <- H_num + H_den
      projector_error <- max(abs(H_test %*% H_test - H_test))
      symmetry_error <- max(abs(H_test - t(H_test)))
      if (!is.finite(projector_error) || !is.finite(symmetry_error) ||
          max(projector_error, symmetry_error) >
            max(100 * .Machine$double.eps, sqrt(.Machine$double.eps))) next
      G_test <- H_test %*% object$design$G %*% H_test
      G_test <- (G_test + t(G_test)) / 2
      ev <- eigen(G_test, symmetric = TRUE, only.values = TRUE)$values
      scale <- max(abs(ev))
      if (!is.finite(scale) || scale <= 0) next
      lambda <- ev[abs(ev) > 1e-8 * scale]
      if (!length(lambda)) next
      out[i] <- .monte_carlo_p(lambda, df_num[i], df_den[i], f_obs[i], n_mc)
    }
  }
  out
}

# Normalise any engine `$table` into a PRIMER-style results table.
.permanova_table <- function(object) {
  fit <- object$fit
  tab <- fit$table
  pick <- function(names) {
    for (nm in names) if (!is.null(tab[[nm]])) return(tab[[nm]])
    rep(NA_real_, nrow(tab))
  }
  source <- as.character(.or(tab[["term"]], .or(tab[["contrast"]], tab[["Source"]])))
  pseudo_f <- pick(c("pseudo_f", "statistic"))
  p_value <- pick(c("p_perm", "p_value", "p_sign", "p_closed", "p_adjusted"))
  result <- data.frame(
    Source = source, df = pick("df"), SS = pick("ss"), MS = pick("ms"),
    `Pseudo-F` = pseudo_f, `P(perm)` = p_value,
    check.names = FALSE, stringsAsFactors = FALSE
  )
  # A grand Total over all observations is only meaningful when the reported
  # sources partition every observation degree of freedom. Aggregated-unit
  # engines (random slope, unequal random intercept, unbalanced nested/hybrid)
  # test in a reduced unit space, so their sources do not sum to `n - 1`.
  partitions_full_data <-
    isTRUE(sum(result$df, na.rm = TRUE) == object$design$n - 1L)
  result$`P(MC)` <- .permanova_mc_column(
    object, result$Source, result$df, pick("denominator_df"), result$`Pseudo-F`,
    partitions_full_data
  )
  used <- .permanova_used(fit)
  result$`Unique perms` <- ifelse(is.finite(result$`Pseudo-F`), used, NA_real_)
  if (partitions_full_data && !any(result$Source == "Total")) {
    total <- data.frame(
      Source = "Total", df = object$design$n - 1L,
      SS = sum(diag(object$design$G)), MS = NA_real_, `Pseudo-F` = NA_real_,
      `P(perm)` = NA_real_, `P(MC)` = NA_real_, `Unique perms` = NA_real_,
      check.names = FALSE
    )
    result <- rbind(result, total)
  }
  result
}

.fmt <- function(x, digits) {
  out <- formatC(x, format = "f", digits = digits)
  out[is.na(x)] <- ""
  out
}

.print_permanova_table <- function(tab) {
  display <- data.frame(
    Source = tab$Source,
    df = ifelse(is.na(tab$df), "", formatC(tab$df, format = "d")),
    SS = .fmt(tab$SS, 3), MS = .fmt(tab$MS, 3),
    `Pseudo-F` = .fmt(tab$`Pseudo-F`, 4),
    `P(perm)` = .fmt(tab$`P(perm)`, 4),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  # P(MC) sits beside P(perm) only when the design supports the asymptotic draw.
  if (!is.null(tab$`P(MC)`) && any(is.finite(tab$`P(MC)`))) {
    display$`P(MC)` <- .fmt(tab$`P(MC)`, 4)
  }
  display$`Unique perms` <- ifelse(is.na(tab$`Unique perms`), "",
                                   formatC(tab$`Unique perms`, format = "d"))
  print(display, row.names = FALSE, right = FALSE)
}

.print_permanova_factors <- function(design) {
  rows <- list()
  for (term in design$fixed_effects) {
    value <- design$data[[term]]
    levels <- if (is.factor(value)) nlevels(droplevels(value)) else NA_integer_
    rows[[length(rows) + 1L]] <- data.frame(
      Name = term, Type = "Fixed", Levels = levels, stringsAsFactors = FALSE
    )
  }
  for (spec in design$random_terms) {
    grouping <- design$groups[[spec$label]]
    rows[[length(rows) + 1L]] <- data.frame(
      Name = spec$label, Type = paste0("Random, ", spec$origin),
      Levels = if (is.null(grouping)) NA_integer_ else nlevels(droplevels(grouping)),
      stringsAsFactors = FALSE
    )
  }
  cat("Factors\n")
  print(do.call(rbind, rows), row.names = FALSE, right = FALSE)
}

#' @export
print.permanova <- function(x, ...) {
  cat("Distance-based PERMANOVA\n")
  cat("Engine: ", x$engine,
      if (!is.na(x$type)) paste0("  (Type ", x$type, " SS)") else "", "\n\n", sep = "")
  .print_permanova_factors(x$design)
  cat("\nPERMANOVA table of results\n")
  .print_permanova_table(.permanova_table(x))
  invisible(x)
}

# EMS coefficient, pseudo-F construction, and components of variation for a
# fixed-effect crossed/factorial model tested over the residual (PRIMER Fig
# 1.44 layout). The coefficient for a term's own component is
# `tr(H_term Z Z') / df`, the same trace-derived EMS coefficient the crossed
# random engine uses; the component estimate is `(MS_term - MS_res) / coef`.
.fixed_component_blocks <- function(object) {
  fit <- object$fit
  design <- object$design
  tab <- fit$table
  residual <- which(tab$term == "Residual")
  tests <- fit$tests
  if (!length(residual) || is.null(tests)) return(NULL)
  ms_res <- tab$ms[[residual]]
  df_res <- tab$df[[residual]]

  ems <- fratio <- components <- list()
  for (term in names(tests)) {
    spec <- tests[[term]]
    projector <- spec$numerator$projector
    df <- spec$numerator$df
    ms <- spec$numerator$ms
    variables <- strsplit(term, ":", fixed = TRUE)[[1L]]
    cells <- interaction(lapply(design$data[variables], factor), drop = TRUE)
    indicator <- stats::model.matrix(~ 0 + cells)
    coefficient <- trace_product(projector, tcrossprod(indicator)) / df
    estimate <- if (abs(coefficient) > 1e-10) (ms - ms_res) / coefficient else NA_real_
    ems[[term]] <- data.frame(
      Source = term,
      EMS = sprintf("1*V(Res) + %s*S(%s)",
                    formatC(coefficient, format = "f", digits = 3), term),
      stringsAsFactors = FALSE
    )
    fratio[[term]] <- data.frame(
      Source = term, Numerator = paste0("1*", term), Denominator = "1*Res",
      Num.df = df, Den.df = df_res, check.names = FALSE, stringsAsFactors = FALSE
    )
    components[[term]] <- data.frame(
      Source = paste0("S(", term, ")"), Estimate = estimate,
      Sq.root = sign(estimate) * sqrt(abs(estimate)),
      check.names = FALSE, stringsAsFactors = FALSE
    )
  }
  list(
    ems = rbind(do.call(rbind, ems),
                data.frame(Source = "Res", EMS = "1*V(Res)", stringsAsFactors = FALSE)),
    fratio = do.call(rbind, fratio),
    components = rbind(
      do.call(rbind, components),
      data.frame(Source = "V(Res)", Estimate = ms_res, Sq.root = sqrt(ms_res),
                 check.names = FALSE)
    )
  )
}

#' @export
summary.permanova <- function(object, ...) {
  print(object)

  fixed_blocks <- if (identical(object$engine, "fixed effects")) {
    .fixed_component_blocks(object)
  } else NULL
  if (!is.null(fixed_blocks)) {
    cat("\nDetails of the expected mean squares (EMS) for the model\n")
    print(fixed_blocks$ems, row.names = FALSE, right = FALSE)
    cat("\nConstruction of Pseudo-F ratio(s) from mean squares\n")
    print(fixed_blocks$fratio, row.names = FALSE, right = FALSE)
    cat("\nEstimates of components of variation\n")
    display <- fixed_blocks$components
    display$Estimate <- formatC(display$Estimate, format = "f", digits = 3)
    display$Sq.root <- formatC(display$Sq.root, format = "f", digits = 3)
    print(display, row.names = FALSE, right = FALSE)
    return(invisible(.permanova_table(object)))
  }

  tests <- .permanova_tests(object$fit)
  denominators <- lapply(tests, function(test) test$denominator$expression)
  denominators <- denominators[!vapply(denominators, is.null, logical(1))]
  if (length(denominators)) {
    cat("\nDenominators\n")
    print(data.frame(
      Source = names(denominators),
      Denominator = unlist(denominators, use.names = FALSE),
      stringsAsFactors = FALSE
    ), row.names = FALSE, right = FALSE)
  }

  components <- .or(object$fit$variance_components$unrestricted,
                    object$fit$covariance$unrestricted)
  # PRIMER's single "components of variation" block combines fixed effect sizes
  # (S(term)) and random variance components (V(term)); the prefix carries the
  # distinction. Covariate-adjusted nested fits use that combined layout.
  fixed_components <- if (identical(object$engine, "balanced nested with covariate")) {
    .nested_covariate_fixed_components(object)
  } else NULL

  if (!is.null(fixed_components) && !is.null(components) && !is.matrix(components)) {
    random_order <- c(setdiff(names(components), "residual"), "residual")
    component_labels <- object$fit$component_labels
    display <- vapply(random_order, function(name) {
      if (identical(name, "residual")) "Res"
      else if (!is.null(component_labels) && name %in% names(component_labels)) component_labels[[name]]
      else name
    }, character(1))
    random_block <- data.frame(
      Source = paste0("V(", display, ")"),
      Estimate = unname(components[random_order]),
      Sq.root = sign(components[random_order]) * sqrt(abs(components[random_order])),
      check.names = FALSE, stringsAsFactors = FALSE
    )
    combined <- rbind(fixed_components, random_block)
    cat("\nEstimates of components of variation\n")
    combined$Estimate <- formatC(combined$Estimate, format = "f", digits = 4)
    combined$Sq.root <- formatC(combined$Sq.root, format = "f", digits = 4)
    print(combined, row.names = FALSE, right = FALSE)
    cat("  S(): fixed-effect size (may be negative); V(): random variance component\n")
  } else if (!is.null(components)) {
    cat("\nEstimates of components of variation\n")
    if (is.matrix(components)) {
      print(round(components, 5))
    } else {
      print(data.frame(
        Source = names(components), Estimate = unname(components),
        `Sq.root` = sign(components) * sqrt(abs(components)),
        check.names = FALSE, stringsAsFactors = FALSE
      ), row.names = FALSE, right = FALSE)
    }
  }
  invisible(.permanova_table(object))
}
