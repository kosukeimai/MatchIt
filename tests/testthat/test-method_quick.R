# Tests for method = "quick" (generalized full matching).
# Structural checks plus snapshot pins recording present behavior. The matching is
# delegated to quickmatch, so the pinned values are a joint property of MatchIt and
# the installed quickmatch; a snapshot failure after a quickmatch upgrade is expected
# rather than a MatchIt regression.

skip_if_not_installed("quickmatch")

data("lalonde", package = "MatchIt")

f <- treat ~ age + educ + race + married + nodegree + re74 + re75
fm <- treat ~ age + educ + re74 + re75

lalonde_sw <- seq(0.5, 2, length.out = nrow(lalonde))

test_that("baseline", {
  m <- matchit(f, data = lalonde, method = "quick")
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)
  expect_matchit_snapshot(m)
})

test_that("estimand='ATC'", {
  m <- matchit(f, data = lalonde, method = "quick", estimand = "ATC")
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)
  expect_matchit_snapshot(m)
})

test_that("estimand='ATE'", {
  m <- matchit(f, data = lalonde, method = "quick", estimand = "ATE")
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)
  expect_matchit_snapshot(m)
})

test_that("the estimand changes the weights", {
  ms <- lapply(c("ATT", "ATC", "ATE"), function(e) {
    matchit(f, data = lalonde, method = "quick", estimand = e)
  })

  expect_not_equal(ms[[1L]]$weights, ms[[2L]]$weights)
  expect_not_equal(ms[[1L]]$weights, ms[[3L]]$weights)
})

test_that("distance caliper", {
  m <- matchit(f, data = lalonde, method = "quick",
               caliper = 0.2, std.caliper = FALSE)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)
  expect_matchit_snapshot(m)
})

test_that("mahvars", {
  m <- matchit(f, data = lalonde, method = "quick",
               mahvars = ~ age + educ + re74 + re75)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)
  expect_matchit_snapshot(m)
})

test_that("full Mahalanobis", {
  m <- matchit(fm, data = lalonde, method = "quick", distance = "mahalanobis")
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)
  expect_matchit_snapshot(m)
})

test_that("exact", {
  m <- matchit(f, data = lalonde, method = "quick", exact = ~ race)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)

  #No subclass may span two exact strata
  keep <- !is.na(m$subclass)
  spans <- tapply(as.character(lalonde$race[keep]), m$subclass[keep],
                  function(x) length(unique(x)))
  expect_true(all(spans == 1L))

  expect_matchit_snapshot(m)
})

test_that("discard", {
  m <- matchit(f, data = lalonde, method = "quick", discard = "both")
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)
  expect_true(any(is.na(m$subclass)))
  expect_matchit_snapshot(m)
})

test_that("s.weights", {
  m <- matchit(f, data = lalonde, method = "quick", s.weights = lalonde_sw)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)
  expect_matchit_snapshot(m)
})

test_that("every subclass contains both treatment groups", {
  m <- matchit(f, data = lalonde, method = "quick")

  keep <- !is.na(m$subclass)
  tab <- table(m$subclass[keep], m$treat[keep])

  expect_true(all(tab[, "0"] > 0L))
  expect_true(all(tab[, "1"] > 0L))
})

test_that("include.obj returns the quickmatch object", {
  m <- matchit(f, data = lalonde, method = "quick", include.obj = TRUE)
  expect_s3_class(m$obj, "qm_matching")
})

test_that("matching improves balance", {
  m <- matchit(f, data = lalonde, method = "quick")
  expect_balance_improved(m)
})

test_that("covariate calipers are an error", {
  expect_err(matchit(f, data = lalonde, method = "quick",
                     caliper = c(age = 2), std.caliper = FALSE),
             'With `method = "quick"`, calipers cannot be placed on covariates.')
})

test_that("mahvars plus a caliper is an error", {
  expect_err(matchit(f, data = lalonde, method = "quick",
                     mahvars = ~ age, caliper = 0.2),
             "a caliper can only be used when")
})

test_that("unused arguments warn and are ignored", {
  expect_wrn(
    matchit(f, data = lalonde, method = "quick", ratio = 2),
    'The argument `ratio` is not used with `method = "quick"` and will be ignored.'
  )

  expect_wrn(
    matchit(f, data = lalonde, method = "quick", antiexact = ~ married),
    'The argument `antiexact` is not used with `method = "quick"` and will be ignored.'
  )
})

test_that("missing values in covariates are an error", {
  lalonde_na <- inject_missingness(lalonde, "educ")

  expect_err(matchit(f, data = lalonde_na, method = "quick"),
             "Missing and non-finite values are not allowed in the covariates")
})

test_that("no unexpected conditions in the baseline call", {
  expect_no_unexpected_warning(matchit(f, data = lalonde, method = "quick"))
})
