# Crossed mixed inference for incomplete, unequally replicated designs.
#
# `permanova_crossed()` requires a complete factorial with equal replication:
# every cell present, every cell the same size. Real crossed designs often are
# neither. A field experiment run at several sites and revisited over time
# typically stops at different times at different sites, leaving whole
# factor-by-factor cells empty, and loses the odd replicate along the way.
#
# The consequence of an empty cell is not just imbalance. A term whose cells are
# partly missing cannot carry its nominal degrees of freedom: with s sites, t
# times and only c of the s*t cells occupied, the site-by-time interaction spans
# c - s - t + 1 dimensions rather than (s-1)(t-1). The design is rank deficient,
# which the fully-crossed engine refuses and the unbalanced fixed engine refuses
# separately. This engine accepts it and reports the estimable dimension.
#
# Three pieces make that work:
#
#   * Sum-to-zero coding. Sums of squares here are marginal -- each term is
#     what the full model explains beyond the model without it. Under indicator
#     coding an interaction's columns span its own margins, so dropping a main
#     effect would remove nothing and every main effect would take zero df.
#     Sum-to-zero contrasts separate the margins from the interaction, which is
#     also the convention PRIMER states it uses for mixed terms.
#
#   * Rank, not arithmetic. Every degree of freedom is a difference of
#     projection ranks rather than a product of level counts, so an incomplete
#     term reports what it can actually estimate.
#
#   * Denominators from expected mean squares. With random factors crossed into
#     the design, a fixed term is generally not tested over the residual but
#     over whichever stratum carries the same random components. Under imbalance
#     no single stratum need match, and the denominator becomes a combination
#     with Satterthwaite degrees of freedom. `.ems_f_construction()` does that
#     matching, and is the same routine validated against PRIMER elsewhere.

# Model matrix for a set of terms under sum-to-zero coding, with the column
# block belonging to each term recorded.
.crossed_incomplete_matrix <- function(labels, data) {
  variables <- unique(unlist(strsplit(labels, ":", fixed = TRUE)))
  missing <- setdiff(variables, names(data))
  if (length(missing)) {
    stop("Unknown design factors: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  frame <- data[variables]
  # Numeric columns stay numeric: a covariate carries one degree of freedom as a
  # slope, and turning it into a factor would spend one per distinct value.
  # Character and logical columns become factors; existing factors are dropped
  # to their observed levels so empty cells cannot inflate the nominal count.
  categorical <- character()
  for (name in variables) {
    if (is.numeric(frame[[name]])) next
    if (!is.factor(frame[[name]])) frame[[name]] <- factor(frame[[name]])
    frame[[name]] <- droplevels(frame[[name]])
    categorical <- c(categorical, name)
  }
  coding <- stats::setNames(rep(list("contr.sum"), length(categorical)), categorical)
  formula <- stats::reformulate(labels)
  structure_terms <- stats::terms(formula, data = frame)
  X <- stats::model.matrix(structure_terms, frame,
                           contrasts.arg = if (length(coding)) coding else NULL)
  assign <- attr(X, "assign")
  ordered <- attr(structure_terms, "term.labels")
  blocks <- lapply(seq_along(ordered), function(k) which(assign == k))
  # `terms()` orders an interaction's components by where its variables first
  # appear in the formula, so a term asked for as `chamber:time` can come back
  # as `time:chamber`. Key the blocks on the sorted variable set so a term is
  # found by what it contains rather than by how it was written.
  names(blocks) <- vapply(ordered, .canonical_term, character(1))
  list(X = X, blocks = blocks, intercept = which(assign == 0L),
       order = vapply(ordered, .canonical_term, character(1)))
}

# Projector for every term, adjusted for all terms that do not contain it.
#
# The tempting construction -- drop the term's columns from the full model and
# see what the rank loses -- fails here. A term that contains this one spans it:
# with the site-by-time interaction still in the model, removing the site
# columns removes nothing the rest cannot supply, and site reports zero
# estimable degrees of freedom. That is not a property of the data; it is the
# wrong comparison. In a complete design the higher-order contrasts are
# orthogonal to the margin and the question does not arise, which is why the
# balanced engine can drop columns safely.
#
# So each term is compared against the model holding every term that does NOT
# contain it, with the terms that do contain it excluded from both sides. Term
# `u` contains term `t` when t's variables are a subset of u's. The comparison
# is then always well posed, and its rank difference is the term's estimable
# dimension whether or not cells are missing.
.crossed_incomplete_projectors <- function(matrix_parts, labels,
                                           hypothesis = "II", tol = 1e-10) {
  X <- matrix_parts$X
  variables <- lapply(stats::setNames(labels, labels), function(label) {
    strsplit(label, ":", fixed = TRUE)[[1L]]
  })
  # Sequential order is the model matrix's own: main effects, then two-way
  # terms, and so on, which is the order a PERMANOVA table is normally printed
  # in.
  sequence <- matrix_parts$order
  canonical <- vapply(labels, .canonical_term, character(1))
  lapply(stats::setNames(labels, labels), function(label) {
    outside <- switch(
      hypothesis,
      # Everything else in the model, including the terms that contain this one.
      "III" = setdiff(labels, label),
      # Everything that does not contain this term.
      "II" = labels[!vapply(labels, function(other) {
        all(variables[[label]] %in% variables[[other]])
      }, logical(1))],
      # Everything earlier in the sequence.
      "I" = labels[canonical %in%
                     sequence[seq_len(match(.canonical_term(label), sequence) - 1L)]],
      stop("`hypothesis` must be \"I\", \"II\" or \"III\".", call. = FALSE)
    )
    reduced_columns <- c(matrix_parts$intercept, unlist(
      matrix_parts$blocks[vapply(outside, .canonical_term, character(1))],
      use.names = FALSE))
    term_columns <- c(reduced_columns,
                      matrix_parts$blocks[[.canonical_term(label)]])
    reduced <- projection_matrix(X[, reduced_columns, drop = FALSE], tol = tol)
    including <- projection_matrix(X[, term_columns, drop = FALSE], tol = tol)
    list(H = including$H - reduced$H, H_reduced = reduced$H,
         df = including$rank - reduced$rank,
         nominal_df = length(matrix_parts$blocks[[.canonical_term(label)]]))
  })
}

# The restricted parameterisation, applied to a fixed term's own expected mean
# square.
#
# `.restricted_ems()` handles the strata; this handles the numerator. Under
# fixed effects summing to zero, the expected mean square of a fixed term A
# carries the random component A-by-random and nothing else: a component
# carrying an additional fixed factor C averages away over C's levels, so
# E(MS_bird) contains sigma^2(bird x time) but not sigma^2(site x bird x time).
#
# A random component is therefore kept only when the fixed factors it involves
# are exactly the term's own. The residual is always carried.
# Which random components appear in a fixed term's expected mean square.
#
# A component belongs there when the term sits inside it and the component
# brings no *other* fixed factor along. The second half has to be judged by
# crossing, not by the label: `chamber(treatment) x time` carries the word
# treatment only because chamber is nested in treatment, and the factor
# genuinely crossing it is time. Reading the label literally would say the
# component belongs to `treatment x time` alone, and would drop it from
# E(MS_time) -- leaving time tested over the residual instead of over
# chamber(treatment) x time.
#
# `crossing` gives, for each component, the fixed factors that actually cross
# it, so the rule is: keep the component when the term's variables sit inside
# the component's, and everything crossing the component is part of the term.
.restrict_fixed_ems <- function(numerator_ems, label, crossing) {
  own <- strsplit(label, ":", fixed = TRUE)[[1L]]
  for (component in names(numerator_ems)) {
    if (identical(component, "residual")) next
    parts <- strsplit(component, ":", fixed = TRUE)[[1L]]
    keep <- all(own %in% parts) && all(crossing[[component]] %in% own)
    if (!keep) numerator_ems[[component]] <- 0
  }
  numerator_ems
}

#' Crossed mixed PERMANOVA for incomplete, unequally replicated designs
#'
#' Fixed and random factors fully crossed, where cells may be empty and
#' replication may be unequal. Sums of squares are marginal: each term is tested
#' for what it adds to the model containing every other term. Denominators come
#' from the design's expected mean squares, so a fixed term crossed with a random
#' factor is tested over the stratum carrying the same random components, and
#' falls back to a Satterthwaite combination when imbalance leaves no single
#' stratum matching.
#'
#' Degrees of freedom are estimable dimensions, not level-count products. A term
#' whose cells are partly empty reports the rank it actually spans, and
#' `nominal_df` in the returned table records what a complete design would have
#' given, so the loss is visible rather than silent.
#'
#' Expected mean squares follow the restricted parameterisation -- fixed effects
#' summing to zero over mixed terms -- so a random stratum carries no component
#' from a random-by-fixed interaction sitting inside it.
#'
#' @param design A `permanova_design` from [compile_permanova_design()] with no
#'   `/`-generated nesting. Random terms are optional: with none, every term is
#'   tested over the residual, which is how an incomplete *fixed* design is
#'   handled.
#' @param hypothesis Sums of squares. `"II"` (default) adjusts each term for
#'   every term that does not contain it. `"I"` is sequential, in the model
#'   matrix's term order. `"III"` adjusts each term for every other term in the
#'   model; on an incomplete design a term containing this one can span it, so
#'   Type III may report fewer estimable degrees of freedom than Type II -- and
#'   in a saturated random model can drive a main effect to zero. Type II is the
#'   default because it stays well posed however incomplete the design is.
#' @param permutations Number of Freedman-Lane permutations, or `"all"`.
#' @param seed Optional random seed.
#' @param keep_permutations Retain the permutation statistics.
#' @param max_exact Largest permutation space to enumerate exactly.
#' @param tol Numerical tolerance for ranks.
#' @return An object of class `permanova_crossed_incomplete`.
#' @seealso [permanova_crossed()] for complete, equally replicated designs.
#' @export
permanova_crossed_incomplete <- function(design, hypothesis = c("II", "I", "III"),
                                         permutations = 999L,
                                         seed = NULL, keep_permutations = FALSE,
                                         max_exact = 100000L, tol = 1e-10) {
  hypothesis <- match.arg(hypothesis)
  if (!inherits(design, "permanova_design")) {
    stop("`design` must be produced by `compile_permanova_design()`.",
         call. = FALSE)
  }
  if (is.null(design$G)) {
    stop("The design object must include `distance`.", call. = FALSE)
  }
  # A random term is not required. With none, every stratum collapses to the
  # residual and each term is tested over it -- an incomplete fixed design,
  # which the Type III fixed engine refuses because the empty cells alias its
  # columns. The rank-based construction here does not mind.
  origins <- if (length(design$random_terms)) {
    vapply(design$random_terms, `[[`, character(1), "origin")
  } else character()
  if (any(origins == "nesting")) {
    stop("This engine covers crossed designs; `/`-generated nesting is handled ",
         "by the nested engines.", call. = FALSE)
  }

  G <- design$G
  n <- nrow(G)
  fixed_labels <- design$fixed_effects
  random_labels <- vapply(design$random_terms, `[[`, character(1), "label")
  if (!length(fixed_labels)) {
    stop("Crossed mixed inference requires at least one fixed term.",
         call. = FALSE)
  }
  labels <- c(fixed_labels, random_labels)

  parts <- .crossed_incomplete_matrix(labels, design$data)
  present <- names(parts$blocks)
  wanted <- vapply(labels, .canonical_term, character(1))
  if (!all(wanted %in% present)) {
    stop("Terms absent from the model matrix: ",
         paste(labels[!wanted %in% present], collapse = ", "), call. = FALSE)
  }
  projectors <- .crossed_incomplete_projectors(parts, labels,
                                              hypothesis = hypothesis, tol = tol)

  full <- projection_matrix(parts$X, tol = tol)
  H_residual <- diag(n) - full$H
  residual_df <- n - full$rank
  if (residual_df <= 0L) {
    stop("The model saturates the data; no residual degrees of freedom remain.",
         call. = FALSE)
  }
  residual_ms <- trace_hgh(H_residual, G) / residual_df

  # A random effect's covariance is Z Z' as declared. It is not orthogonalised
  # against the terms containing it: an inner component genuinely appears in an
  # outer stratum's expected mean square, and the restricted parameterisation --
  # not a change of basis -- is what removes it where PRIMER removes it.
  covariance <- c(
    stats::setNames(lapply(random_labels, function(label) {
      tcrossprod(design$Z[[label]])
    }), random_labels),
    list(residual = diag(n))
  )
  if (!length(random_labels)) covariance <- list(residual = diag(n))
  components <- names(covariance)
  ems_of <- function(H, df) {
    vapply(covariance, function(V) trace_product(H, V) / df, numeric(1))
  }

  strata_labels <- c(random_labels, "Residual")
  strata_H <- c(lapply(projectors[random_labels], `[[`, "H"), list(H_residual))
  names(strata_H) <- strata_labels
  strata_df <- c(vapply(projectors[random_labels], `[[`, numeric(1), "df"),
                 Residual = residual_df)
  usable <- strata_df > 0
  if (!all(usable)) {
    dropped <- strata_labels[!usable]
    warning("Random stratum with no estimable dimension dropped: ",
            paste(dropped, collapse = ", "), call. = FALSE)
  }
  strata_labels <- strata_labels[usable]
  strata_H <- strata_H[usable]
  strata_df <- strata_df[usable]
  strata_ms <- vapply(strata_labels, function(label) {
    trace_hgh(strata_H[[label]], G) / strata_df[[label]]
  }, numeric(1))

  candidate <- do.call(rbind, lapply(strata_labels, function(label) {
    ems_of(strata_H[[label]], strata_df[[label]])
  }))
  rownames(candidate) <- strata_labels
  # Random terms that also involve a fixed factor are random-by-fixed
  # interactions, which the restricted parameterisation removes from the
  # expected mean square of any stratum sitting inside them.
  fixed_variables <- unique(unlist(strsplit(fixed_labels, ":", fixed = TRUE)))
  # A random term counts as a random-by-fixed interaction only where a fixed
  # factor genuinely crosses it. Carrying a fixed factor's name is not enough:
  # nesting is written the same way. In `treatment:chamber` the treatment is
  # fixed by the chamber -- each chamber sits in one treatment -- so the term is
  # a nested random effect, not an interaction. In `treatment:chamber:time` the
  # time varies within every chamber, so that one is. The test is whether the
  # fixed factor is determined by the term's other variables.
  determined_by <- function(target, others) {
    others <- others[!vapply(design$data[others], is.numeric, logical(1))]
    if (!length(others) || is.numeric(design$data[[target]])) return(FALSE)
    key <- interaction(design$data[others], drop = TRUE)
    all(tapply(as.character(design$data[[target]]), key,
               function(z) length(unique(z))) == 1L)
  }
  crossing_fixed <- lapply(stats::setNames(nm = random_labels), function(label) {
    variables <- strsplit(label, ":", fixed = TRUE)[[1L]]
    candidates <- intersect(variables, fixed_variables)
    candidates[vapply(candidates, function(f) {
      !determined_by(f, setdiff(variables, f))
    }, logical(1))]
  })
  interaction_terms <- random_labels[lengths(crossing_fixed) > 0L]
  candidate <- .restricted_ems(candidate, interaction_terms)

  if (!is.null(seed)) set.seed(seed)
  plan <- .permutation_orders(n, permutations, max_exact)

  rows <- list()
  permutation_statistics <- list()
  for (label in labels) {
    spec <- projectors[[label]]
    row <- data.frame(
      term = label, df = spec$df, nominal_df = spec$nominal_df,
      ss = trace_hgh(spec$H, G), ms = NA_real_, pseudo_f = NA_real_,
      denominator = NA_character_, denominator_df = NA_real_,
      p_perm = NA_real_, stringsAsFactors = FALSE
    )
    if (spec$df <= 0L) {
      rows[[length(rows) + 1L]] <- row
      next
    }
    row$ms <- row$ss / spec$df
    # A random term's numerator is one of the strata, so it takes the restricted
    # row already computed rather than a fresh unrestricted one: under the
    # restricted parameterisation E(MS_time) carries no component from the
    # random-by-fixed interactions sitting inside it, which is what leaves it
    # tested over the residual. A fixed term keeps every component it carries --
    # E(MS_A) does contain the A-by-random interaction -- so it is computed
    # directly.
    numerator_ems <- if (label %in% strata_labels) {
      candidate[label, ]
    } else {
      .restrict_fixed_ems(ems_of(spec$H, spec$df), label, crossing_fixed)
    }
    # The null being tested is that this term contributes nothing, so the
    # denominator must have the numerator's expected mean square WITHOUT the
    # term's own contribution. For a fixed term that is automatic: its effect
    # enters the EMS as a quadratic form, not as one of the variance components
    # being matched. For a random term it is not -- the term is itself one of
    # the candidate strata, and matching it against itself would return the
    # term over the term and a pseudo-F of exactly one. So the component under
    # test is zeroed out of the numerator and its own stratum withheld from the
    # candidates.
    available <- seq_along(strata_labels)
    if (label %in% strata_labels) {
      numerator_ems[[label]] <- 0
      available <- setdiff(available, match(label, strata_labels))
    }
    built <- .ems_f_construction(numerator_ems, row$ms, spec$df,
                                 candidate[available, , drop = FALSE],
                                 strata_ms[available], strata_df[available],
                                 tol = tol)
    if (is.null(built)) {
      rows[[length(rows) + 1L]] <- row
      next
    }
    # Restore full-length weights so the denominator label and the permutation
    # both index the complete stratum list.
    coefficients <- numeric(length(strata_labels))
    coefficients[available] <- built$coefficients
    built$coefficients <- coefficients
    row$pseudo_f <- built$statistic
    row$denominator_df <- built$denominator_df
    row$denominator <- paste(sprintf(
      "%.5g*%s", built$coefficients[built$coefficients > tol],
      strata_labels[built$coefficients > tol]), collapse = " + ")
    if (nrow(plan$orders)) {
      statistics <- .freedman_lane_quasi_statistics(
        G, spec$H_reduced, spec$H, spec$df, strata_H, strata_df,
        built$coefficients, plan$orders,
        label = paste0("Permuting for ", label),
        numerator_residual = built$numerator_residual,
        residual_index = match("Residual", strata_labels), tol = tol
      )
      inference <- .permutation_inference(statistics, built$statistic,
                                          plan$exact, tol)
      row$p_perm <- inference$p_value
      if (keep_permutations) permutation_statistics[[label]] <- statistics
    }
    rows[[length(rows) + 1L]] <- row
  }

  table <- do.call(rbind, rows)
  table <- rbind(table, data.frame(
    term = "Residual", df = residual_df, nominal_df = NA_integer_,
    ss = trace_hgh(H_residual, G), ms = residual_ms, pseudo_f = NA_real_,
    denominator = NA_character_, denominator_df = NA_real_, p_perm = NA_real_,
    stringsAsFactors = FALSE
  ), data.frame(
    term = "Total", df = n - 1L, nominal_df = NA_integer_,
    ss = sum(diag(G)), ms = NA_real_, pseudo_f = NA_real_,
    denominator = NA_character_, denominator_df = NA_real_, p_perm = NA_real_,
    stringsAsFactors = FALSE
  ))
  rownames(table) <- NULL

  incomplete <- table$term %in% labels &
    !is.na(table$nominal_df) & table$df < table$nominal_df
  out <- list(
    call = match.call(), table = table, design = design,
    fixed_effects = fixed_labels, random_terms = random_labels,
    ems = candidate, strata_ms = strata_ms, strata_df = strata_df,
    hypothesis = hypothesis,
    permutations = if (is.character(permutations)) permutations else
      as.integer(permutations),
    exact = plan$exact, seed = seed,
    rank_deficient = table$term[incomplete],
    permutation_statistics = if (keep_permutations) permutation_statistics else NULL,
    method = paste("marginal SS; expected-mean-square denominators;",
                   "Freedman-Lane permutation of residuals")
  )
  class(out) <- "permanova_crossed_incomplete"
  out
}

#' @export
print.permanova_crossed_incomplete <- function(x, ...) {
  cat("Crossed mixed PERMANOVA (incomplete / unequally replicated)\n")
  cat("Method: Type ", x$hypothesis, " ", x$method, "\n", sep = "")
  cat("Random terms:", paste(x$random_terms, collapse = ", "), "\n\n")
  show <- x$table
  show$ss <- signif(show$ss, 5)
  show$ms <- signif(show$ms, 5)
  show$pseudo_f <- round(show$pseudo_f, 4)
  show$denominator_df <- round(show$denominator_df, 2)
  print(show, row.names = FALSE)
  if (length(x$rank_deficient)) {
    cat("\nEstimable df below nominal (empty cells):",
        paste(x$rank_deficient, collapse = ", "), "\n")
  }
  invisible(x)
}
