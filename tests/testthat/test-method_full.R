test_that("calipers work, positive and negative", {
  set.seed(1234)
  n <- 1e3
  p <- runif(n, 0, .4)
  x <- matrix(runif(n * 4), nrow = n)
  g <- sample(1:5, n, TRUE)
  a <- rbinom(n, 1, p)
  dis <- as.logical(rbinom(n, 1, .1))
  d <- data.frame(p, x, a, g, dis)

  #Positive calipers
  m <- matchit(a ~ X1 + X2 + X3 + X4, data = d, distance = d$p, method = "full",
               caliper = .001, std.caliper = FALSE)

  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)

  m <- matchit(a ~ X1 + X2 + X3 + X4, data = d, distance = d$p, method = "full",
               caliper = c(X1 = .001),
               std.caliper = FALSE)

  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)

  m <- matchit(a ~ X1 + X2 + X3 + X4, data = d, distance = d$p, method = "full",
               caliper = c(.02, X1 = .01),
               std.caliper = FALSE)

  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)

  #Negative calipers
  m <- matchit(a ~ X1 + X2 + X3 + X4, data = d, distance = d$p, method = "full",
               caliper = -.3,
               std.caliper = FALSE)

  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)

  m <- matchit(a ~ X1 + X2 + X3 + X4, data = d, distance = d$p, method = "full",
               caliper = c(X1 = -.5),
               std.caliper = FALSE)

  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)

  m <- matchit(a ~ X1 + X2 + X3 + X4, data = d, distance = d$p, method = "full",
               caliper = c(-.3, X1 = -.5),
               std.caliper = FALSE)

  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)

  m <- matchit(a ~ X1 + X2 + X3 + X4, data = d, distance = d$p, method = "full",
               caliper = c(-.3, X1 = .001),
               std.caliper = FALSE)

  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)

  m <- matchit(a ~ X1 + X2 + X3 + X4, data = d, distance = d$p, method = "full",
               caliper = c(.002, X1 = -.5),
               std.caliper = FALSE)

  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)
})
# ===========================================================================
# Full matching: option surface on lalonde, with pinned results.
#
# The matching is delegated to optmatch, so the pinned values are a joint property
# of MatchIt and the installed optmatch; a snapshot failure after an optmatch
# upgrade is expected rather than a MatchIt regression. All specifications below
# are deterministic.
# ===========================================================================

skip_if_not_installed("optmatch")

data("lalonde", package = "MatchIt")

f_full <- treat ~ age + educ + race + married + nodegree + re74 + re75
fm_full <- treat ~ age + educ + re74 + re75

lalonde_sw <- seq(0.5, 2, length.out = nrow(lalonde))

#Full matching retains every unit unless something removes it, and every stratum must
#contain at least one unit from each treatment group.
expect_full_strata_valid <- function(m) {
  keep <- !is.na(m$subclass)

  tab <- table(m$subclass[keep], m$treat[keep])
  expect_true(all(tab[, "0"] > 0L))
  expect_true(all(tab[, "1"] > 0L))

  invisible(m)
}

test_that("full: baseline retains every unit", {
  m <- matchit(f_full, data = lalonde, method = "full")
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)
  expect_full_strata_valid(m)
  expect_false(anyNA(m$subclass))
  expect_true(all(m$weights > 0))
  expect_matchit_snapshot(m)
})

test_that("full: estimand", {
  for (e in c("ATT", "ATC", "ATE")) {
    m <- matchit(f_full, data = lalonde, method = "full", estimand = e)
    expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                        expect_subclass = TRUE)
    expect_full_strata_valid(m)
    expect_matchit_snapshot(m)
  }
})

test_that("full: the estimand changes the weights", {
  ms <- lapply(c("ATT", "ATC", "ATE"), function(e) {
    matchit(f_full, data = lalonde, method = "full", estimand = e)
  })

  expect_not_equal(ms[[1L]]$weights, ms[[2L]]$weights)
  expect_not_equal(ms[[1L]]$weights, ms[[3L]]$weights)
})

test_that("full: min.controls", {
  m <- matchit(f_full, data = lalonde, method = "full", min.controls = 1)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)

  #Every stratum must now hold at least one control per treated unit
  keep <- !is.na(m$subclass)
  tab <- table(m$subclass[keep], m$treat[keep])
  expect_true(all(tab[, "0"] >= tab[, "1"]))

  expect_matchit_snapshot(m)
})

test_that("full: max.controls", {
  m <- matchit(f_full, data = lalonde, method = "full", max.controls = 4)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)

  #No stratum may hold more than four controls per treated unit
  keep <- !is.na(m$subclass)
  tab <- table(m$subclass[keep], m$treat[keep])
  expect_true(all(tab[, "0"] <= 4L * tab[, "1"]))

  expect_matchit_snapshot(m)
})

test_that("full: min.controls and max.controls together", {
  m <- matchit(f_full, data = lalonde, method = "full",
               min.controls = 0.5, max.controls = 4)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)

  keep <- !is.na(m$subclass)
  tab <- table(m$subclass[keep], m$treat[keep])
  expect_true(all(tab[, "0"] >= 0.5 * tab[, "1"]))
  expect_true(all(tab[, "0"] <= 4L * tab[, "1"]))

  expect_matchit_snapshot(m)
})

test_that("full: tighter restrictions produce more strata", {
  #Restricting the ratio forces the optimizer to split the sample more finely.
  n_sub <- function(...) {
    nlevels(matchit(f_full, data = lalonde, method = "full", ...)$subclass)
  }

  expect_gt(n_sub(min.controls = 1), n_sub())
  expect_gt(n_sub(max.controls = 4), n_sub())
})

test_that("full: mean.controls drops units to hit the average", {
  m <- matchit(f_full, data = lalonde, method = "full", mean.controls = 2)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)
  expect_true(any(is.na(m$subclass)))
  expect_matchit_snapshot(m)
})

test_that("full: omit.fraction drops units", {
  m <- matchit(f_full, data = lalonde, method = "full", omit.fraction = 0.2)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)
  expect_true(any(is.na(m$subclass)))
  expect_lt(sum(m$weights > 0), nrow(lalonde))
  expect_matchit_snapshot(m)
})

test_that("full: mean.controls and omit.fraction cannot both be given", {
  expect_error(matchit(f_full, data = lalonde, method = "full",
                       mean.controls = 2, omit.fraction = 0.2),
               "cannot both be specified")
})

test_that("full: distance caliper", {
  m <- matchit(f_full, data = lalonde, method = "full",
               caliper = 0.2, std.caliper = FALSE)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)
  expect_matchit_snapshot(m)
})

test_that("full: a caliper cannot be combined with min.controls", {
  expect_err(matchit(f_full, data = lalonde, method = "full",
                     caliper = 0.2, min.controls = 1),
             'Calipers cannot be used with `method = "full"` when `min.controls` is specified.')
})

test_that("full: exact", {
  m <- matchit(f_full, data = lalonde, method = "full", exact = ~ race)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)

  #No stratum may span two levels of the exact matching variable
  keep <- !is.na(m$subclass)
  spans <- tapply(as.character(lalonde$race[keep]), m$subclass[keep],
                  function(x) length(unique(x)))
  expect_true(all(spans == 1L))

  expect_matchit_snapshot(m)
})

test_that("full: antiexact", {
  m <- matchit(f_full, data = lalonde, method = "full", antiexact = ~ married)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)
  expect_matchit_snapshot(m)
})

test_that("full: mahvars", {
  m <- matchit(f_full, data = lalonde, method = "full", mahvars = ~ age + educ)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)
  expect_matchit_snapshot(m)
})

test_that("full: full Mahalanobis", {
  m <- matchit(fm_full, data = lalonde, method = "full", distance = "mahalanobis")
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)
  expect_matchit_snapshot(m)
})

test_that("full: a supplied distance matrix", {
  lalonde_sub <- lalonde[c(1:50, 186:235), ]
  d <- scaled_euclidean_dist(fm_full, data = lalonde_sub)

  m <- matchit(fm_full, data = lalonde_sub, method = "full", distance = d)
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)
  expect_matchit_snapshot(m)
})

test_that("full: discard", {
  m <- matchit(f_full, data = lalonde, method = "full", discard = "both")
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)
  expect_true(any(is.na(m$subclass)))
  expect_matchit_snapshot(m)
})

test_that("full: s.weights", {
  m <- matchit(f_full, data = lalonde, method = "full", s.weights = lalonde_sw)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)
  expect_matchit_snapshot(m)
})

test_that("full: tol is passed through to optmatch", {
  m <- matchit(f_full, data = lalonde, method = "full", tol = 1e-7)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)

  #A tighter tolerance changes the solution the solver settles on
  expect_not_equal(m$subclass,
                   matchit(f_full, data = lalonde, method = "full")$subclass)

  expect_matchit_snapshot(m)
})

test_that("full: include.obj returns the optmatch object", {
  m <- matchit(f_full, data = lalonde, method = "full", include.obj = TRUE)
  expect_s3_class(m$obj, "optmatch")
})

test_that("full: full matching is no worse than optimal pair matching on total distance", {
  #Full matching minimizes the total within-stratum distance over a strictly larger
  #set of allowed structures than 1:1 optimal matching, so it cannot do worse.
  mf <- matchit(f_full, data = lalonde, method = "full")
  mo <- matchit(f_full, data = lalonde, method = "optimal", tol = 1e-6)

  mean_within <- function(m) {
    keep <- !is.na(m$subclass)
    mean(vapply(levels(droplevels(m$subclass[keep])), function(s) {
      i <- which(keep & m$subclass == s)
      mean(abs(outer(m$distance[i[m$treat[i] == 1L]],
                     m$distance[i[m$treat[i] == 0L]], "-")))
    }, numeric(1L)))
  }

  expect_lte(mean_within(mf), mean_within(mo))
})

test_that("full: matching improves balance", {
  expect_balance_improved(matchit(f_full, data = lalonde, method = "full"))
})

test_that("full: unused arguments warn and are ignored", {
  m0 <- matchit(f_full, data = lalonde, method = "full")

  for (arg in c("replace", "ratio", "m.order")) {
    args <- list(f_full, data = lalonde, method = "full")
    args[[arg]] <- switch(arg, replace = TRUE, ratio = 2, m.order = "data")

    expect_wrn(do.call(matchit, args),
               sprintf('The argument `%s` is not used with `method = "full"` and will be ignored.',
                       arg))

    expect_identical(suppressWarnings(do.call(matchit, args))$subclass, m0$subclass)
  }
})

test_that("full: mahvars with a full-distance `distance` is an error", {
  expect_err(matchit(fm_full, data = lalonde, method = "full",
                     distance = "mahalanobis", mahvars = ~ age),
             "cannot be used with")
})

test_that("full: missing values in covariates are an error", {
  lalonde_na <- inject_missingness(lalonde, "educ")

  expect_err(matchit(f_full, data = lalonde_na, method = "full"),
             "Missing and non-finite values are not allowed in the covariates")
})

test_that("full: no unexpected conditions in the baseline call", {
  expect_no_unexpected_warning(matchit(f_full, data = lalonde, method = "full"))
})
