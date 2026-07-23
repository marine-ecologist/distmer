# Build an EMS table (sources x variance components) for a 1-D response by the
# trace method, so the structured constructor can be checked end to end.
.ems_table <- function(d, y, term_matrices, component_matrices, leading) {
  n <- nrow(d)
  G <- gower_center(dist(y))
  proj <- function(M) {
    if (!ncol(M)) return(matrix(0, n, n))
    q <- qr(M, tol = 1e-10)
    Q <- qr.Q(q)[, seq_len(q$rank), drop = FALSE]
    tcrossprod(Q)
  }
  cumulative <- matrix(1, n, 1)
  previous <- proj(cumulative)
  H <- list()
  for (name in names(term_matrices)) {
    cumulative <- cbind(cumulative, term_matrices[[name]])
    current <- proj(cumulative)
    H[[name]] <- current - previous
    previous <- current
  }
  H[["Residual"]] <- diag(n) - previous
  df <- vapply(H, function(h) as.integer(round(sum(diag(h)))), integer(1))
  ms <- vapply(H, function(h) sum(h * t(G)), numeric(1)) / df
  ems <- t(vapply(names(H), function(name) {
    vapply(component_matrices, function(V) sum(H[[name]] * t(V)) / df[name], numeric(1))
  }, numeric(length(component_matrices))))
  colnames(ems) <- names(component_matrices)
  rownames(ems) <- names(H)
  list(ems = ems, mean_squares = ms, df = df)
}

test_that("clean balanced nested design yields exact single-source denominators", {
  set.seed(1)
  d <- expand.grid(rep = 1:4, Area = 1:2, Site = 1:2, Location = 1:3)
  d$Location <- factor(d$Location)
  d$Site <- factor(paste0(d$Location, d$Site))
  d$Area <- factor(paste0(d$Site, d$Area))
  d$y <- 2 * (d$Location == "2") + rnorm(nlevels(d$Site), 0, 1.5)[d$Site] +
    rnorm(nlevels(d$Area), 0, 1)[d$Area] + rnorm(nrow(d))
  mm <- function(f) model.matrix(f, d)[, -1, drop = FALSE]
  Zs <- model.matrix(~ 0 + Site, d); Za <- model.matrix(~ 0 + Area, d)
  tab <- .ems_table(
    d, d$y,
    list(Location = mm(~ Location), Site = mm(~ Site), Area = mm(~ Area)),
    list(residual = diag(nrow(d)), site = tcrossprod(Zs), area = tcrossprod(Za)),
    leading = c(Site = "site", Area = "area", Residual = "residual")
  )
  leading <- c(Site = "site", Area = "area", Residual = "residual")

  loc <- structured_denominator(tab$ems, tab$mean_squares, tab$df, leading, "Location")
  sit <- structured_denominator(tab$ems, tab$mean_squares, tab$df, leading, "Site")
  are <- structured_denominator(tab$ems, tab$mean_squares, tab$df, leading, "Area")

  expect_equal(names(loc$weights), "Site")       # Location tests over Site(Location)
  expect_equal(unname(loc$weights), 1, tolerance = 1e-8)
  expect_equal(names(sit$weights), "Area")       # Site tests over Area(Site)
  expect_equal(names(are$weights), "Residual")   # Area tests over the residual
  expect_equal(loc$df, tab$df[["Site"]], tolerance = 1e-8)
})

test_that("random-slope denominator df agrees with lmerTest on a centred covariate", {
  skip_if_not_installed("lmerTest")
  set.seed(3)
  d <- data.frame(Site = factor(rep(1:12, each = 6)))
  d$Volume <- rnorm(nrow(d), 0, 2)
  d$Vw <- decompose_covariate(d, "Volume", by = "Site")$Volume_within
  intercept <- rnorm(12, 0, 1.5)[d$Site]; slope <- rnorm(12, 0, 0.6)[d$Site]
  d$y <- 1 + 0.5 * d$Vw + intercept + slope * d$Vw + rnorm(nrow(d), 0, 1)

  reference <- anova(lmerTest::lmer(
    y ~ Vw + (1 | Site) + (0 + Vw | Site), d, REML = TRUE
  ))
  ref_ddf <- reference[["DenDF"]][1]

  Zs <- model.matrix(~ 0 + Site, d); Dv <- diag(d$Vw)
  tab <- .ems_table(
    d, d$y,
    list(Vo = matrix(d$Vw, ncol = 1), Site = model.matrix(~ Site, d)[, -1, drop = FALSE],
         VoSite = d$Vw * model.matrix(~ Site, d)[, -1, drop = FALSE]),
    list(residual = diag(nrow(d)), site = tcrossprod(Zs),
         voSite = Dv %*% tcrossprod(Zs) %*% Dv),
    leading = c(Site = "site", VoSite = "voSite", Residual = "residual")
  )
  leading <- c(Site = "site", VoSite = "voSite", Residual = "residual")
  vo <- structured_denominator(tab$ems, tab$mean_squares, tab$df, leading, "Vo")
  # the centred slope is tested primarily over the random-slope term; a small
  # residual correction remains because the fixed- and slope-term VoSite EMS
  # coefficients differ slightly
  expect_true("VoSite" %in% names(vo$weights))
  expect_equal(names(vo$weights)[which.max(abs(vo$weights))], "VoSite")
  # Satterthwaite df roughly agrees with lmerTest (method-of-moments vs REML)
  expect_lt(abs(vo$df - ref_ddf), 3)
})

test_that("structured_denominator validates its inputs", {
  ems <- rbind(Site = c(residual = 1, site = 4), Residual = c(residual = 1, site = 0))
  expect_error(
    structured_denominator(ems, c(Site = 3), c(Site = 5, Residual = 40),
                           c(Residual = "residual", Site = "site"), "Site"),
    "named for every source"
  )
})
