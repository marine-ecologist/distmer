.unwrap_parens <- function(x) {
  while (is.call(x) && identical(x[[1L]], as.name("("))) x <- x[[2L]]
  x
}

.is_random_term <- function(x) {
  x <- .unwrap_parens(x)
  is.call(x) && identical(x[[1L]], as.name("|"))
}

.collect_random_terms <- function(x) {
  if (.is_random_term(x)) return(list(.unwrap_parens(x)))
  if (!is.call(x)) return(list())
  unlist(lapply(as.list(x)[-1L], .collect_random_terms), recursive = FALSE)
}

.drop_random_terms <- function(x) {
  if (.is_random_term(x)) return(NULL)
  if (!is.call(x)) return(x)

  op <- as.character(x[[1L]])
  if (op == "+") {
    left <- .drop_random_terms(x[[2L]])
    right <- .drop_random_terms(x[[3L]])
    if (is.null(left)) return(right)
    if (is.null(right)) return(left)
    return(call("+", left, right))
  }

  if (length(.collect_random_terms(x))) {
    stop("Random-effect terms must be top-level additive terms.", call. = FALSE)
  }
  x
}

.group_atomic_variables <- function(x) {
  x <- .unwrap_parens(x)
  if (is.name(x)) return(as.character(x))
  if (is.call(x) && identical(x[[1L]], as.name(":"))) {
    return(c(.group_atomic_variables(x[[2L]]), .group_atomic_variables(x[[3L]])))
  }
  stop("Grouping interactions may contain only variable names and `:`.", call. = FALSE)
}

.parse_group_expression <- function(x) {
  x <- .unwrap_parens(x)
  if (is.name(x)) {
    variable <- as.character(x)
    return(list(list(
      label = variable, variables = variable, origin = "main",
      parent = NA_character_, expression = variable
    )))
  }
  if (is.call(x) && identical(x[[1L]], as.name(":"))) {
    variables <- unique(.group_atomic_variables(x))
    return(list(list(
      label = paste(variables, collapse = ":"), variables = variables,
      origin = "interaction", parent = NA_character_,
      expression = paste(deparse(x), collapse = "")
    )))
  }
  if (is.call(x) && identical(x[[1L]], as.name("/"))) {
    parent_specs <- .parse_group_expression(x[[2L]])
    child_variables <- unique(.group_atomic_variables(x[[3L]]))
    parent <- parent_specs[[length(parent_specs)]]
    variables <- unique(c(parent$variables, child_variables))
    child <- list(
      label = paste(variables, collapse = ":"), variables = variables,
      origin = "nesting", parent = parent$label,
      expression = paste(deparse(x), collapse = "")
    )
    return(c(parent_specs, list(child)))
  }
  stop("Grouping expressions may contain only variable names, `:`, and `/`.", call. = FALSE)
}

.parse_random_slope_left <- function(x) {
  x <- .unwrap_parens(x)
  if (!is.call(x) || !identical(x[[1L]], as.name("+"))) {
    stop(
      "Random-slope terms must use the explicit form `(1 + slope | group)`.",
      call. = FALSE
    )
  }
  terms <- as.list(x)[-1L]
  intercept <- vapply(terms, function(term) {
    identical(term, 1) || identical(term, 1L)
  }, logical(1))
  slopes <- terms[!intercept]
  if (sum(intercept) != 1L || length(slopes) != 1L ||
      !is.name(.unwrap_parens(slopes[[1L]]))) {
    stop(
      "Random-slope terms must contain one intercept and one named slope: `(1 + slope | group)`.",
      call. = FALSE
    )
  }
  as.character(.unwrap_parens(slopes[[1L]]))
}

.parse_random_slope_spec <- function(left, right) {
  slope <- .parse_random_slope_left(left)
  right <- .unwrap_parens(right)
  if (!is.name(right)) {
    stop("Phase 27 random slopes require one named grouping factor.",
         call. = FALSE)
  }
  group <- as.character(right)
  list(
    label = paste0("1+", slope, "|", group), variables = group,
    origin = "slope", parent = NA_character_,
    expression = paste0("1 + ", slope, " | ", group),
    slope = slope, correlated = TRUE, group = group
  )
}

.deduplicate_group_specs <- function(specs) {
  if (!length(specs)) return(specs)
  labels <- vapply(specs, `[[`, character(1), "label")
  retained <- list()
  for (label in unique(labels)) {
    matching <- specs[labels == label]
    signatures <- vapply(matching, function(x) {
      slope <- if (is.null(x$slope)) "" else x$slope
      paste(x$origin, x$parent, paste(x$variables, collapse = ":"), slope,
            sep = "|")
    }, character(1))
    if (length(unique(signatures)) > 1L) {
      stop(
        "Grouping term `", label,
        "` is declared with conflicting nesting/interaction provenance.",
        call. = FALSE
      )
    }
    retained[[length(retained) + 1L]] <- matching[[1L]]
  }
  retained
}

.random_slope_matrix <- function(spec, data) {
  missing <- setdiff(c(spec$group, spec$slope), names(data))
  if (length(missing)) {
    stop("Unknown random-slope variable(s): ",
         paste(missing, collapse = ", "), ".", call. = FALSE)
  }
  group <- droplevels(factor(data[[spec$group]]))
  group_matrix <- stats::model.matrix(~ 0 + group)
  slope_frame <- data.frame(.slope = data[[spec$slope]])
  slope_matrix <- stats::model.matrix(~ .slope, slope_frame)[, -1L, drop = FALSE]
  if (!ncol(slope_matrix)) {
    stop("The random-slope variable must contain estimable variation.",
         call. = FALSE)
  }
  slope_columns <- do.call(cbind, lapply(seq_len(ncol(slope_matrix)), function(j) {
    group_matrix * as.numeric(slope_matrix[, j])
  }))
  matrix <- cbind(group_matrix, slope_columns)
  colnames(matrix) <- c(
    paste0(spec$group, "[", levels(group), "]:intercept"),
    unlist(lapply(colnames(slope_matrix), function(name) {
      paste0(spec$group, "[", levels(group), "]:", name)
    }))
  )
  matrix
}

.indicator_matrix <- function(path, data) {
  variables <- strsplit(path, ":", fixed = TRUE)[[1L]]
  missing <- setdiff(variables, names(data))
  if (length(missing)) {
    stop("Unknown grouping variable(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  values <- lapply(data[variables], factor)
  group <- if (length(values) == 1L) values[[1L]] else interaction(values, drop = TRUE)
  matrix <- stats::model.matrix(~ 0 + group)
  colnames(matrix) <- paste0(path, "[", levels(group), "]")
  matrix
}

.fixed_term_projections <- function(X, fixed_terms, tol) {
  assignment <- attr(X, "assign")
  labels <- attr(fixed_terms, "term.labels")
  full <- projection_matrix(X, tol = tol)

  setNames(lapply(seq_along(labels), function(k) {
    reduced_X <- X[, assignment != k, drop = FALSE]
    reduced <- projection_matrix(reduced_X, tol = tol)
    list(
      term = labels[[k]],
      columns = which(assignment == k),
      H_reduced = reduced$H,
      H_term = full$H - reduced$H,
      df = full$rank - reduced$rank,
      estimable = full$rank > reduced$rank
    )
  }), labels)
}

.crossing_table <- function(specs) {
  if (length(specs) < 2L) {
    return(data.frame(
      first = character(), second = character(), crossed = logical(),
      relation = character()
    ))
  }
  pairs <- utils::combn(seq_along(specs), 2L, simplify = FALSE)
  do.call(rbind, lapply(pairs, function(pair) {
    first <- specs[[pair[[1L]]]]
    second <- specs[[pair[[2L]]]]
    nested <- (first$origin == "nesting" && identical(first$parent, second$label)) ||
      (second$origin == "nesting" && identical(second$parent, first$label))
    component <- !nested && (all(first$variables %in% second$variables) ||
      all(second$variables %in% first$variables))
    crossed <- !nested && !component &&
      !length(intersect(first$variables, second$variables))
    relation <- if (nested) "nested" else if (component) "component" else if (crossed) {
      "crossed"
    } else "overlapping"
    data.frame(
      first = first$label, second = second$label, crossed = crossed,
      relation = relation
    )
  }))
}

#' Compile a canonical distance-model design object
#'
#' Parses fixed effects, random-intercept terms, and the narrow correlated
#' random-slope syntax used by Phase 27 without performing any
#' inference. The resulting object is intended to be the sole model-structure
#' input to later EMS and permutation engines.
#'
#' @param formula An `lme4`-style formula. Random effects support
#'   `(1 | group)`, nesting with `/`, explicit grouping interactions with `:`,
#'   and one correlated `(1 + slope | group)` declaration.
#' @param data A data frame with one row per observation.
#' @param distance Optional distance object used for row-count and spectral
#'   diagnostics.
#' @param tol Numerical rank and eigenvalue tolerance.
#' @return An object of class `permanova_design`.
#' @export
compile_permanova_design <- function(formula, data, distance = NULL, tol = 1e-10) {
  if (!inherits(formula, "formula")) stop("`formula` must be a formula.", call. = FALSE)
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  if (nrow(data) < 2L) stop("`data` must contain at least two observations.", call. = FALSE)

  rhs <- formula[[length(formula)]]
  random_calls <- .collect_random_terms(rhs)
  fixed_rhs <- .drop_random_terms(rhs)
  if (is.null(fixed_rhs)) fixed_rhs <- 1

  random_specs <- list()
  for (term in random_calls) {
    intercept <- .unwrap_parens(term[[2L]])
    if (identical(intercept, 1) || identical(intercept, 1L)) {
      random_specs <- c(random_specs, .parse_group_expression(term[[3L]]))
    } else {
      random_specs <- c(
        random_specs, list(.parse_random_slope_spec(intercept, term[[3L]]))
      )
    }
  }
  random_specs <- .deduplicate_group_specs(random_specs)
  random_paths <- vapply(random_specs, `[[`, character(1), "label")

  fixed_formula <- as.formula(call("~", fixed_rhs), env = environment(formula))
  fixed_terms <- terms(fixed_formula, data = data)
  fixed_frame <- model.frame(fixed_terms, data = data, na.action = na.fail,
                             drop.unused.levels = TRUE)
  X <- model.matrix(fixed_terms, fixed_frame)
  fixed_projection <- projection_matrix(X, tol = tol)
  term_projections <- .fixed_term_projections(X, fixed_terms, tol)

  random_matrices <- setNames(lapply(random_specs, function(spec) {
    if (identical(spec$origin, "slope")) {
      .random_slope_matrix(spec, data)
    } else {
      .indicator_matrix(spec$label, data)
    }
  }), random_paths)
  random_projections <- lapply(random_matrices, projection_matrix, tol = tol)
  fixed_aliases <- if (fixed_projection$rank < ncol(X)) {
    pivot <- qr(X, tol = tol)$pivot
    colnames(X)[pivot[(fixed_projection$rank + 1L):ncol(X)]]
  } else {
    character()
  }
  nested_specs <- Filter(function(x) x$origin == "nesting", random_specs)
  containment <- if (length(nested_specs)) {
    do.call(rbind, lapply(nested_specs, function(x) data.frame(
      parent = x$parent, child = x$label, stringsAsFactors = FALSE
    )))
  } else data.frame(parent = character(), child = character())

  crossing <- .crossing_table(random_specs)
  grouping <- setNames(lapply(random_specs, function(spec) {
    variables <- spec$variables
    values <- lapply(data[variables], factor)
    if (length(values) == 1L) values[[1L]] else interaction(values, drop = TRUE)
  }), random_paths)
  random_term_table <- if (length(random_specs)) {
    data.frame(
      term = random_paths,
      origin = vapply(random_specs, `[[`, character(1), "origin"),
      parent = vapply(random_specs, function(x) x$parent, character(1)),
      variables = vapply(random_specs, function(x) paste(x$variables, collapse = ":"), character(1)),
      expression = vapply(random_specs, `[[`, character(1), "expression"),
      slope = vapply(random_specs, function(x) {
        if (is.null(x$slope)) NA_character_ else x$slope
      }, character(1)),
      correlated = vapply(random_specs, function(x) {
        isTRUE(x$correlated)
      }, logical(1)),
      stringsAsFactors = FALSE
    )
  } else data.frame(
    term = character(), origin = character(), parent = character(),
    variables = character(), expression = character(), slope = character(),
    correlated = logical()
  )
  group_sizes <- lapply(grouping, table)
  balance <- list(
    random_group_sizes = group_sizes,
    random_balanced = vapply(group_sizes, function(x) length(unique(as.integer(x))) <= 1L, logical(1))
  )

  spectral <- NULL
  D <- NULL
  G <- NULL
  if (!is.null(distance)) {
    D <- as_distance_matrix(distance)
    if (nrow(D) != nrow(data)) {
      stop("`data` rows must match the distance-matrix observation order.", call. = FALSE)
    }
    spectral <- distance_diagnostics(D, tol = tol)
    G <- spectral$G
    spectral$G <- NULL
  }

  response <- if (length(formula) == 3L) paste(deparse(formula[[2L]]), collapse = "") else NULL
  warnings <- character()
  if (length(fixed_aliases)) warnings <- c(warnings, "Fixed-effect model matrix is rank deficient.")
  if (!is.null(spectral) && !spectral$euclidean) {
    warnings <- c(warnings, "Dissimilarity is non-Euclidean; negative eigenvalues were retained.")
  }

  out <- list(
    call = match.call(), formula = formula, fixed_formula = fixed_formula,
    response = response, data = data,
    fixed_effects = attr(fixed_terms, "term.labels"), random_effects = random_paths,
    random_terms = random_specs, random_term_table = random_term_table,
    nesting = containment, crossing = crossing,
    interaction_graph = list(
      fixed = attr(fixed_terms, "factors"),
      random = random_term_table[random_term_table$origin == "interaction", , drop = FALSE]
    ),
    containment_graph = containment, X = X, Z = random_matrices, groups = grouping,
    projections = list(fixed = fixed_projection, terms = term_projections, random = random_projections),
    estimable_functions = lapply(term_projections, function(x) x$H_term),
    estimability = list(
      fixed_rank = fixed_projection$rank, fixed_columns = ncol(X),
      aliased = fixed_aliases,
      terms = vapply(term_projections, `[[`, logical(1), "estimable")
    ),
    exchangeability_units = list(
      observation = seq_len(nrow(data)), random = grouping,
      candidates = c("observation", random_paths)
    ),
    permutation_strata = grouping, balance = balance,
    distance = D, G = G, spectral = spectral, n = nrow(data),
    warnings = warnings, status = "design compiled; inference not specified"
  )
  class(out) <- "permanova_design"
  out
}

#' @export
print.permanova_design <- function(x, ...) {
  cat("Canonical PERMANOVA design\n")
  cat("Observations:", x$n, " Fixed rank:", x$estimability$fixed_rank, "\n")
  cat("Fixed:", if (length(x$fixed_effects)) paste(x$fixed_effects, collapse = ", ") else "(intercept only)", "\n")
  cat("Random:", if (length(x$random_effects)) paste(x$random_effects, collapse = ", ") else "(none)", "\n")
  if (!is.null(x$spectral)) cat("Euclidean distance:", x$spectral$euclidean, "\n")
  if (length(x$warnings)) cat("Warnings:", paste(x$warnings, collapse = " "), "\n")
  invisible(x)
}
