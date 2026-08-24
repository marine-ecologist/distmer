# Multi-level unbalanced nested inference. The engine aggregates by cascading
# equal weights up the chain, so every check below is pinned against an
# independently hand-cascaded reference rather than against the engine itself.

# Build a hierarchy from an explicit shape spec so topology is controlled:
# `shapes[[region]][[site]]` is the vector of per-quadrat replicate counts.
.multi_fixture <- function(shapes, treatment, seed = 9) {
  rows <- list()
  for (r in seq_along(shapes)) {
    for (s in seq_along(shapes[[r]])) {
      counts <- shapes[[r]][[s]]
      for (q in seq_along(counts)) {
        rows[[length(rows) + 1L]] <- data.frame(
          reg = r, site = s, quad = q, rep = seq_len(counts[[q]]),
          tr = treatment[[r]], stringsAsFactors = FALSE
        )
      }
    }
  }
  d <- do.call(rbind, rows)
  for (v in c("reg", "site", "quad", "tr")) d[[v]] <- factor(d[[v]])
  set.seed(seed)
  d$y <- 1.2 * (d$tr == "t") + stats::rnorm(nlevels(d$reg))[d$reg] +
    stats::rnorm(nrow(d), 0, 0.6)
  d
}

# Two shapes, each shared by several regions, so topology strata span the
# treatment contrast while replication stays unequal at two levels.
.shape_a <- list(c(3, 3), c(3, 3))
.shape_b <- list(c(2), c(2, 2, 2))

# Hand cascade, computed with tapply only. Unit labelling is taken from the
# compiled design so the comparison is order-aligned; the averaging itself is
# independent of the engine.
.hand_cascade <- function(d, design) {
  quad <- droplevels(factor(design$groups[["reg:site:quad"]]))
  site <- droplevels(factor(design$groups[["reg:site"]]))
  region <- droplevels(factor(design$groups[["reg"]]))
  quad_mean <- tapply(d$y, quad, mean)
  site_mean <- tapply(quad_mean, site[match(names(quad_mean), quad)], mean)
  region_mean <- tapply(site_mean, region[match(names(site_mean), site)], mean)
  list(quad = quad_mean, site = site_mean, region = region_mean,
       quad_factor = quad, site_factor = site, region_factor = region)
}

test_that("cascaded membership reproduces hand aggregation at every level", {
  d <- .multi_fixture(list(.shape_a, .shape_b, .shape_a, .shape_a, .shape_b, .shape_a),
                      c("c", "c", "c", "t", "t", "t"))
  design <- compile_permanova_design(~ tr + (1 | reg / site / quad), d, dist(d$y))
  chain <- .nested_chain(design)
  expect_equal(chain$depth, 3L)

  hand <- .hand_cascade(d, design)
  expected <- list(hand$region, hand$site, hand$quad)
  for (level in seq_len(chain$depth)) {
    cascade <- .cascade_membership(chain, level)
    expect_equal(as.vector(cascade$M %*% d$y),
                 as.vector(expected[[level]]), tolerance = 1e-10)
    # every unit averages its own children with weights summing to one
    expect_equal(unname(rowSums(cascade$M)),
                 rep(1, nrow(cascade$M)), tolerance = 1e-12)
  }
})

test_that("top-level treatment matches an lm on hand-cascaded region centroids", {
  d <- .multi_fixture(list(.shape_a, .shape_b, .shape_a, .shape_a, .shape_b, .shape_a),
                      c("c", "c", "c", "t", "t", "t"))
  fit <- permanova_unbalanced_nested_multi(
    compile_permanova_design(~ tr + (1 | reg / site / quad), d, dist(d$y)),
    permutations = 0
  )
  design <- compile_permanova_design(~ tr + (1 | reg / site / quad), d, dist(d$y))
  hand <- .hand_cascade(d, design)
  treatment <- factor(tapply(as.character(d$tr), hand$region_factor,
                             function(v) v[[1L]]))
  reference <- anova(lm(hand$region ~ treatment))

  expect_equal(fit$tests[["tr"]]$hierarchy_level, "reg")
  expect_equal(fit$table$pseudo_f[[1L]], reference[1L, "F value"], tolerance = 1e-10)
  expect_equal(fit$table$df[[1L]], reference[1L, "Df"])
  expect_equal(fit$table$denominator_df[[1L]], reference[2L, "Df"])
})

test_that("lower-level treatment is tested within its parent level", {
  shape <- list(c(3, 3), c(3, 3, 3), c(3, 3), c(3, 3, 3))
  d <- .multi_fixture(rep(list(shape), 4L), rep("c", 4L), seed = 11)
  # treatment varies across sites within region -> level 2
  d$tr <- factor(ifelse(as.integer(d$site) <= 2L, "c", "t"))
  set.seed(11)
  d$y <- 0.9 * (d$tr == "t") + stats::rnorm(4)[d$reg] + stats::rnorm(nrow(d), 0, 0.5)

  fit <- permanova_unbalanced_nested_multi(
    compile_permanova_design(~ tr + (1 | reg / site / quad), d, dist(d$y)),
    permutations = 0
  )
  design <- compile_permanova_design(~ tr + (1 | reg / site / quad), d, dist(d$y))
  hand <- .hand_cascade(d, design)
  region <- hand$region_factor[match(names(hand$site), hand$site_factor)]
  treatment <- d$tr[match(names(hand$site), hand$site_factor)]
  reference <- anova(lm(hand$site ~ region + treatment))

  expect_equal(fit$tests[["tr"]]$hierarchy_level, "reg:site")
  expect_equal(fit$table$pseudo_f[[1L]],
               reference["treatment", "F value"], tolerance = 1e-10)
  expect_match(fit$tests[["tr"]]$permutation_strata, "reg")
})

test_that("permutation never moves a unit outside its parent", {
  shape <- list(c(3, 3), c(3, 3, 3), c(3, 3), c(3, 3, 3))
  d <- .multi_fixture(rep(list(shape), 4L), rep("c", 4L), seed = 11)
  d$tr <- factor(ifelse(as.integer(d$site) <= 2L, "c", "t"))
  set.seed(11)
  d$y <- stats::rnorm(nrow(d))
  design <- compile_permanova_design(~ tr + (1 | reg / site / quad), d, dist(d$y))
  chain <- .nested_chain(design)
  aggregated <- .cascade_aggregate(design, chain, 2L)
  parent_factor <- factor(vapply(
    split(as.character(chain$groups[[1L]]), aggregated$unit), `[[`, character(1), 1L))
  signature <- .nested_topology_signature(chain, 2L)
  strata <- droplevels(interaction(parent_factor, signature, drop = TRUE, sep = " | "))
  set.seed(1)
  plan <- .stratified_permutation_plan(strata, 200L, 100000L)

  parent <- as.integer(parent_factor)
  for (row in seq_len(nrow(plan$orders))) {
    expect_identical(parent[plan$orders[row, ]], parent)
  }
})

test_that("a four-level chain is supported", {
  rows <- list()
  for (r in 1:6) for (s in 1:2) for (q in 1:2) {
    for (b in seq_len(if (r %% 2L) 2L else 3L)) {
      rows[[length(rows) + 1L]] <- data.frame(
        reg = r, site = s, quad = q, sub = b, rep = seq_len(2L))
    }
  }
  d <- do.call(rbind, rows)
  for (v in c("reg", "site", "quad", "sub")) d[[v]] <- factor(d[[v]])
  d$tr <- factor(ifelse(as.integer(d$reg) <= 3L, "c", "t"))
  set.seed(3)
  d$y <- stats::rnorm(6)[d$reg] + stats::rnorm(nrow(d), 0, 0.5)

  design <- compile_permanova_design(~ tr + (1 | reg / site / quad / sub), d, dist(d$y))
  expect_equal(.nested_chain(design)$depth, 4L)
  fit <- permanova_unbalanced_nested_multi(design, permutations = 99L, seed = 1L)
  expect_equal(fit$tests[["tr"]]$hierarchy_level, "reg")
  expect_true(is.finite(fit$table$pseudo_f[[1L]]))
})

test_that("non-exchangeable topologies and observation-level terms are refused", {
  # every region a different shape -> no two regions exchangeable
  shapes <- list(list(c(2)), list(c(3)), list(c(4)), list(c(5)))
  d <- .multi_fixture(shapes, c("c", "c", "t", "t"))
  expect_error(
    permanova_unbalanced_nested_multi(
      compile_permanova_design(~ tr + (1 | reg / site / quad), d, dist(d$y)),
      permutations = 0),
    "do not span every contrast"
  )

  d2 <- .multi_fixture(list(.shape_a, .shape_b, .shape_a, .shape_a, .shape_b, .shape_a),
                       c("c", "c", "c", "t", "t", "t"))
  # a term varying inside the deepest unit is now an observation-level term,
  # tested over the within-unit residual rather than refused
  d2$tr <- factor(rep_len(c("c", "t"), nrow(d2)))
  observation_fit <- permanova_unbalanced_nested_multi(
    compile_permanova_design(~ tr + (1 | reg / site / quad), d2, dist(d2$y)),
    permutations = 0)
  expect_equal(observation_fit$tests[["tr"]]$kind, "observation")
  quad <- factor(paste(d2$reg, d2$site, d2$quad, sep = "."))
  reference <- anova(lm(d2$y ~ quad + d2$tr))
  expect_equal(observation_fit$table$pseudo_f[[1L]],
               reference["d2$tr", "F value"], tolerance = 1e-10)
})

test_that("permanova() routes depth-3 chains here and leaves parent/child alone", {
  d <- .multi_fixture(list(.shape_a, .shape_b, .shape_a, .shape_a, .shape_b, .shape_a),
                      c("c", "c", "c", "t", "t", "t"))
  deep <- compile_permanova_design(~ tr + (1 | reg / site / quad), d, dist(d$y))
  expect_equal(.permanova_dispatch(deep, NULL)$engine,
               "unbalanced nested (multi-level)")

  # a two-level chain keeps its existing validated engine
  shallow_rows <- do.call(rbind, lapply(1:6, function(p)
    do.call(rbind, lapply(seq_len(if (p %% 2L) 2L else 3L), function(c_)
      data.frame(par = p, ch = c_, rep = seq_len(3L))))))
  shallow_rows$par <- factor(shallow_rows$par)
  shallow_rows$ch <- factor(shallow_rows$ch)
  shallow_rows$tr <- factor(ifelse(as.integer(shallow_rows$par) <= 3L, "c", "t"))
  set.seed(2)
  shallow_rows$y <- stats::rnorm(nrow(shallow_rows))
  shallow <- compile_permanova_design(~ tr + (1 | par / ch), shallow_rows,
                                      dist(shallow_rows$y))
  expect_equal(.permanova_dispatch(shallow, NULL)$engine,
               "unbalanced nested (top unit)")
})

test_that("the two-level and multi-level engines share one estimand", {
  # Both must average child centroids with equal weight, so a parent/child
  # design gives the same pseudo-F from either engine and from an lm on
  # hand-cascaded parent centroids.
  sizes <- list(c(5, 2), c(5, 2), c(5, 2), c(2, 5), c(2, 5), c(2, 5))
  rows <- list()
  for (p in seq_along(sizes)) {
    for (child in seq_along(sizes[[p]])) {
      rows[[length(rows) + 1L]] <- data.frame(
        par = p, ch = child, rep = seq_len(sizes[[p]][[child]]),
        tr = if (p <= 3L) "control" else "impact", stringsAsFactors = FALSE
      )
    }
  }
  d <- do.call(rbind, rows)
  for (v in c("par", "ch", "tr")) d[[v]] <- factor(d[[v]])
  set.seed(8)
  d$y <- stats::rnorm(nlevels(d$par))[d$par] + stats::rnorm(nrow(d), 0, 0.5)
  design <- compile_permanova_design(~ tr + (1 | par / ch), d, dist(d$y))

  two_level <- permanova_unbalanced_nested(design, permutations = 0)
  multi <- permanova_unbalanced_nested_multi(design, permutations = 0)

  child_label <- factor(paste(d$par, d$ch, sep = "."))
  child_mean <- tapply(d$y, child_label, mean)
  parent_mean <- tapply(child_mean, d$par[match(names(child_mean), child_label)], mean)
  treatment <- factor(tapply(as.character(d$tr), d$par, function(v) v[[1L]]))
  reference <- anova(lm(parent_mean ~ treatment))[1L, "F value"]

  expect_equal(two_level$table$pseudo_f[[1L]], reference, tolerance = 1e-10)
  expect_equal(multi$table$pseudo_f[[1L]], reference, tolerance = 1e-10)

  # and it must NOT equal the observation-weighted parent centroid here
  observation_weighted <- anova(lm(
    tapply(d$y, d$par, mean) ~ treatment))[1L, "F value"]
  expect_false(isTRUE(all.equal(reference, observation_weighted)))
})

test_that("multiple fixed terms are each tested at their own hierarchy level", {
  # A, B and A:B live at region level; C at site level; x inside the deepest
  # unit. Every term is pinned to its own independent lm reference.
  set.seed(21)
  big <- c(1, 2, 5, 6)
  rows <- list()
  for (r in 1:8) {
    for (s in 1:2) {
      for (q in seq_len(if (r %in% big) 3L else 2L)) {
        rows[[length(rows) + 1L]] <- data.frame(
          reg = r, site = s, quad = q, rep = seq_len(2L))
      }
    }
  }
  d <- do.call(rbind, rows)
  for (v in c("reg", "site", "quad")) d[[v]] <- factor(d[[v]])
  d$A <- factor(ifelse(as.integer(d$reg) %% 2L == 1L, "a1", "a2"))
  d$B <- factor(ifelse(as.integer(d$reg) <= 4L, "b1", "b2"))
  d$C <- factor(ifelse(as.integer(d$site) == 1L, "c1", "c2"))
  d$x <- stats::rnorm(nrow(d))
  d$y <- 1.0 * (d$A == "a2") + 0.6 * (d$C == "c2") + 0.5 * d$x +
    stats::rnorm(8)[d$reg] + stats::rnorm(nrow(d), 0, 0.5)

  fit <- permanova_unbalanced_nested_multi(
    compile_permanova_design(~ A * B + C + x + (1 | reg / site / quad), d,
                             dist(d$y)),
    permutations = 0
  )
  expect_equal(unname(fit$term_levels[c("A", "B", "C", "x", "A:B")]),
               c(1L, 1L, 2L, NA_integer_, 1L))

  quad <- factor(paste(d$reg, d$site, d$quad, sep = "."))
  site <- factor(paste(d$reg, d$site, sep = "."))
  quad_mean <- tapply(d$y, quad, mean)
  site_mean <- tapply(quad_mean, site[match(names(quad_mean), quad)], mean)
  region_mean <- tapply(site_mean, d$reg[match(names(site_mean), site)], mean)

  region <- anova(lm(region_mean ~
    d$A[match(names(region_mean), d$reg)] * d$B[match(names(region_mean), d$reg)]))
  site_ref <- anova(lm(site_mean ~ d$reg[match(names(site_mean), site)] +
    d$C[match(names(site_mean), site)]))
  observation <- anova(lm(d$y ~ quad + d$x))

  value <- function(term) fit$table$pseudo_f[fit$table$term == term]
  expect_equal(value("A"), region[1L, "F value"], tolerance = 1e-9)
  expect_equal(value("B"), region[2L, "F value"], tolerance = 1e-9)
  expect_equal(value("A:B"), region[3L, "F value"], tolerance = 1e-9)
  expect_equal(value("C"), site_ref[2L, "F value"], tolerance = 1e-9)
  expect_equal(value("x"), observation[2L, "F value"], tolerance = 1e-9)

  # one stratum row per occupied level, plus the residual
  expect_true(all(c("reg stratum", "reg:site stratum", "Residual") %in%
                    fit$table$term))
})

test_that("sum-of-squares type changes the partition when cells are unbalanced", {
  set.seed(31)
  A <- c(rep("a1", 6), rep("a2", 6))
  B <- c("b1", "b1", "b1", "b2", "b2", "b2", "b1", "b1", "b2", "b2", "b2", "b2")
  big <- c(1, 2, 4, 5, 7, 8, 9, 10)
  rows <- list()
  for (r in 1:12) {
    for (s in 1:2) {
      for (q in seq_len(if (r %in% big) 3L else 2L)) {
        rows[[length(rows) + 1L]] <- data.frame(
          reg = r, site = s, quad = q, rep = seq_len(2L))
      }
    }
  }
  d <- do.call(rbind, rows)
  for (v in c("reg", "site", "quad")) d[[v]] <- factor(d[[v]])
  d$A <- factor(A[as.integer(d$reg)])
  d$B <- factor(B[as.integer(d$reg)])
  d$y <- 0.8 * (d$A == "a2") + 0.5 * (d$B == "b2") +
    stats::rnorm(12)[d$reg] + stats::rnorm(nrow(d), 0, 0.5)

  quad <- factor(paste(d$reg, d$site, d$quad, sep = "."))
  site <- factor(paste(d$reg, d$site, sep = "."))
  quad_mean <- tapply(d$y, quad, mean)
  site_mean <- tapply(quad_mean, site[match(names(quad_mean), quad)], mean)
  centroid <- data.frame(
    y = as.vector(tapply(site_mean, d$reg[match(names(site_mean), site)], mean)),
    A = factor(A), B = factor(B)
  )

  old <- options(contrasts = c("contr.sum", "contr.poly"))
  on.exit(options(old), add = TRUE)
  design <- compile_permanova_design(~ A * B + (1 | reg / site / quad), d,
                                     dist(d$y))
  type_one <- permanova_unbalanced_nested_multi(design, permutations = 0,
                                                hypothesis = "I")
  type_three <- permanova_unbalanced_nested_multi(design, permutations = 0,
                                                  hypothesis = "III")
  sequential <- anova(lm(y ~ A * B, centroid))

  expect_equal(type_one$table$pseudo_f[1:3],
               unname(sequential[1:3, "F value"]), tolerance = 1e-9)
  # the partition must actually change for a main effect
  expect_false(isTRUE(all.equal(type_one$table$pseudo_f[[1L]],
                                type_three$table$pseudo_f[[1L]])))
  skip_if_not_installed("car")
  partial <- car::Anova(lm(y ~ A * B, centroid), type = 3)
  expect_equal(type_three$table$pseudo_f[1:3],
               unname(partial[c("A", "B", "A:B"), "F value"]), tolerance = 1e-9)
})
