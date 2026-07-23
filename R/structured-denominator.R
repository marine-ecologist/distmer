#' Structured (Cornfield--Tukey / Satterthwaite) denominator construction
#'
#' Builds a pseudo-F denominator for a mixed-model term by *synthesising* mean
#' squares in containment-lattice order, rather than by the unconstrained
#' minimum-norm match used in [quasi_denominator()]. Given the expected mean
#' squares (EMS) of a set of sources, each random source's variance component,
#' and the observed mean squares, it peels the target EMS one variance component
#' at a time from the top of the lattice down, using the source that directly
#' estimates each component. On a clean balanced design this returns the exact
#' single-source denominator; on a perturbed design (e.g. an unbalanced or
#' covariate-containing mixed model) it returns the structured linear
#' combination with Satterthwaite degrees of freedom, matching the construction
#' used by classical mixed-model software.
#'
#' @param ems Numeric matrix of EMS coefficients, one row per source (named) and
#'   one column per variance component (named). Entry `ems[i, j]` is the
#'   coefficient of component `j` in the expected mean square of source `i`.
#' @param mean_squares Named numeric vector of observed mean squares, one per
#'   source (names matching `rownames(ems)`).
#' @param dfs Named numeric vector of degrees of freedom, one per source.
#' @param leading Named character vector mapping each *random* source to the
#'   variance component it directly estimates (e.g. `c(Residual = "residual",
#'   Site = "site")`); fixed sources are omitted or `NA`.
#' @param target The source to test: the denominator is built to match its EMS
#'   with its own leading component removed. May instead be a named numeric
#'   target EMS vector.
#' @param tol Numerical tolerance for treating coefficients as zero.
#' @return A list with the denominator `weights` (named), its value, and its
#'   Satterthwaite `df`.
#' @examples
#' # A balanced two-level nested design: Site over Residual is exact.
#' ems <- rbind(Site = c(residual = 1, site = 4),
#'              Residual = c(residual = 1, site = 0))
#' structured_denominator(ems, c(Site = 3, Residual = 1), c(Site = 5, Residual = 40),
#'                        leading = c(Residual = "residual", Site = "site"),
#'                        target = "Site")
#' @export
structured_denominator <- function(ems, mean_squares, dfs, leading, target,
                                   tol = 1e-8) {
  ems <- as.matrix(ems)
  sources <- rownames(ems)
  components <- colnames(ems)
  if (is.null(sources) || is.null(components)) {
    stop("`ems` must have source row names and component column names.",
         call. = FALSE)
  }
  if (!setequal(names(mean_squares), sources) || !setequal(names(dfs), sources)) {
    stop("`mean_squares` and `dfs` must be named for every source in `ems`.",
         call. = FALSE)
  }
  mean_squares <- mean_squares[sources]
  dfs <- dfs[sources]
  leading <- leading[!is.na(leading)]
  if (!all(leading %in% components) || anyDuplicated(leading)) {
    stop("`leading` must map sources to distinct variance components.",
         call. = FALSE)
  }
  generator <- setNames(names(leading), leading)

  if (length(target) == 1L && is.character(target)) {
    if (!target %in% sources) stop("`target` is not a source.", call. = FALSE)
    tvec <- ems[target, ]
    if (target %in% names(leading)) tvec[leading[[target]]] <- 0
  } else {
    if (is.null(names(target)) || !setequal(names(target), components)) {
      stop("A numeric `target` must be named for every variance component.",
           call. = FALSE)
    }
    tvec <- target[components]
  }

  # Peel components from the top of the containment lattice downward: a
  # component appearing in fewer sources is more senior and cancelled first.
  seniority <- colSums(abs(ems) > tol)
  order <- names(sort(seniority))

  weights <- setNames(numeric(length(sources)), sources)
  for (component in order) {
    if (abs(tvec[[component]]) <= tol) next
    source <- generator[[component]]
    if (is.null(source)) {
      stop("No random source estimates component `", component,
           "`; the denominator cannot be synthesised.", call. = FALSE)
    }
    weight <- tvec[[component]] / ems[source, component]
    weights[[source]] <- weights[[source]] + weight
    tvec <- tvec - weight * ems[source, ]
  }
  active <- abs(weights) > tol
  contributions <- weights[active] * mean_squares[active]
  value <- sum(contributions)
  satterthwaite <- if (abs(value) <= tol) NA_real_ else {
    value^2 / sum(contributions^2 / dfs[active])
  }
  list(weights = weights[active], value = value, df = satterthwaite)
}
