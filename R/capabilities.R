#' Report the validated inference support matrix
#'
#' Returns one row for every primary inferential engine in the validated support
#' contract. The registry is descriptive: it does not relax the runtime design
#' validators, which remain authoritative for a particular compiled design.
#'
#' @return A data frame containing design envelopes, inferential scales,
#'   denominators, permutation contracts, error-control claims, and exclusions.
#' @export
permanova_capabilities <- function() {
  data.frame(
    engine = c(
      "permanova_oneway", "permanova_fixed", "permanova_mixed",
      "permanova_crossed", "permanova_nested", "permanova_dispersion",
      "permanova_hybrid", "permanova_unbalanced_hybrid",
      "permanova_unbalanced",
      "permanova_unequal_random_intercept",
      "permanova_unbalanced_nested", "permanova_unbalanced_nested_lower",
      "permanova_unbalanced_block",
      "permanova_random_slope",
      "permanova_contrasts", "permanova_marginal_contrasts",
      "permanova_adjusted_contrasts", "permanova_mixed_contrasts",
      "permanova_maxT", "permanova_pairwise_strong",
      "permanova_mixed_maxT"
    ),
    companion = c(
      "", "", "", "", "", "", "", "", "", "", "", "", "", "",
      "", "permanova_marginal_pairwise",
      "permanova_adjusted_pairwise", "permanova_mixed_pairwise",
      "permanova_maxT_pairwise", "",
      "permanova_mixed_maxT_pairwise"
    ),
    design_envelope = c(
      "one categorical fixed factor",
      "fixed sequential model",
      "balanced single random intercept",
      "balanced crossed random intercepts with declared fixed-by-random components",
      "one balanced nested random-intercept hierarchy",
      "one categorical target in a supported location-design structure",
      "balanced treatment crossed with region/site hierarchy",
      "treatment crossed with one unbalanced region/site hierarchy",
      "full-rank unbalanced fixed model",
      "one cluster-level treatment and one unequal random intercept",
      "one top-level treatment in an unbalanced parent/child hierarchy",
      "one child- or observation-level treatment in an unbalanced parent/child hierarchy",
      "incomplete or unequally replicated randomized block",
      "balanced complete two-condition correlated random intercept and slope",
      "one categorical fixed factor, unequal group sizes allowed",
      "full-rank categorical fixed model",
      "full-rank fixed model with factors and ordinary numeric covariates",
      "supported balanced or unbalanced mixed treatment structure",
      "one categorical fixed factor",
      "one categorical fixed factor with independent observations",
      "supported balanced or unbalanced mixed treatment structure"
    ),
    estimand = c(
      "factor centroids", "sequential fixed-term location",
      "equally weighted whole-plot centroids",
      "declared fixed-term location over interaction covariance",
      "fixed-term location at one hierarchy stratum",
      "mean distance to group median or centroid",
      "equally weighted treatment-by-region profiles",
      "equal-site-weight observed treatment-by-region profiles after region adjustment",
      "explicit Type I, II, or III fixed-term hypothesis",
      "equally weighted cluster centroids",
      "equally weighted intact top-unit centroids",
      "equally weighted child centroids or child-adjusted observations",
      "equally weighted observed treatment-by-block cell centroids",
      "population mean complete-subject condition-difference vector",
      "planned or pairwise factor-centroid contrast",
      "equal- or frequency-weighted categorical marginal contrast",
      "contrast at explicit numeric covariate reference values",
      "planned or pairwise experimental-unit contrast",
      "fixed one-factor contrast family",
      "all pairwise complete-distribution equality hypotheses",
      "mixed experimental-unit contrast family"
    ),
    denominator = c(
      "within-group residual", "full-model residual",
      "whole-plot residual EMS", "declared exact or EMS-matched strata",
      "exact hierarchy-level EMS stratum", "design-selected EMS stratum",
      "treatment-by-region interaction",
      "profile-centroid residual after region and treatment",
      "full-model residual",
      "unit-centroid residual", "top-unit-centroid residual",
      "child-centroid residual after parent or observation residual after child",
      "cell-centroid residual after block and treatment",
      "between-subject random-slope difference residual",
      "within-group residual", "full-model residual",
      "full-model residual", "design-selected unit residual or interaction",
      "within-group residual", "within-group residual",
      "design-selected unit residual or interaction"
    ),
    permutation_contract = c(
      "unrestricted raw observations", "term-specific residual permutation",
      "whole-plot residual permutation",
      "canonical within-block cells where supported; otherwise parametric opt-in only",
      "intact hierarchy units within declared parents",
      "full-model residuals under the selected design restrictions",
      "intact treatment profiles within regions",
      "intact profile centroids within region-by-exact-topology strata",
      "term-specific unrestricted Freedman-Lane residuals",
      "whole clusters within exact-size strata",
      "intact top units within topology strata",
      "child centroids within parent-by-size strata or observations within child",
      "cell centroids within block-by-size strata",
      "whole-subject sign flips under central symmetry",
      "contrast-specific unrestricted Freedman-Lane residuals",
      "contrast-specific unrestricted Freedman-Lane residuals",
      "contrast-specific unrestricted Freedman-Lane residuals",
      "contrast-specific unit residuals within design strata",
      "shared raw observations under the complete null",
      "intersection-specific raw labels within connected equality components",
      "shared raw experimental units within design strata under the complete null"
    ),
    error_control = c(
      rep("single-hypothesis permutation p-value", 14L),
      rep("explicit Holm, Bonferroni, BH, BY, or none", 4L),
      "weak FWER under the complete equal-centroid null",
      "strong FWER for pairwise complete-distribution nulls by closure",
      "weak FWER under the complete treatment-null randomization hypothesis"
    ),
    explicit_exclusions = c(
      "mixed or multifactor models",
      "random effects and implicit marginal hypotheses",
      "imbalance, random slopes, crossed or nested structures",
      "undeclared fixed-by-random covariance; general crossed permutations",
      "imbalance or hypotheses spanning hierarchy strata",
      "multiple targets and unsupported design structures",
      "imbalance and omitted treatment-by-hierarchy components",
      "partial site profiles, disconnected incidence, topology-confounded treatment, and unbalanced EMS claims",
      "random effects, aliases, and non-estimable hypotheses",
      "treatment varying within cluster or size-confounded treatment",
      "lower-level treatment or topology-confounded treatment",
      "top-level treatment, disconnected incidence, confounded strata, or multiple nesting trees",
      "disconnected incidence or block-by-size-confounded treatment",
      "more than two conditions, missing or unequal cells, asymmetric difference distributions, and random-slope contrasts",
      "mixed, multifactor, aliased, or invalid contrast definitions",
      "numeric covariates, random effects, and aliased reference grids",
      "implicit covariate values, random effects, aliases, and ill-conditioned fixed models",
      "random slopes, multiple nesting trees, and unsupported mixed structures",
      "strong subset-null control and mixed designs",
      "centroid-only nulls with heterogeneous distributions, planned contrasts, dependent or mixed designs",
      "strong subset-null control, random slopes, and unsupported mixed structures"
    ),
    stringsAsFactors = FALSE
  )
}
