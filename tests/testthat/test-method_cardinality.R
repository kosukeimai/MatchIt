# Tests for method = "cardinality".
# Structural checks plus snapshot pins recording present behavior. The matching is a
# MILP solved by an external solver, so the pinned values are a joint property of
# MatchIt and the installed solver.
#
# Only `solver = "highs"` results are pinned; `"glpk"` is checked for agreement on the
# size of the optimum rather than the identity of the selected units.

skip_if_not_installed("highs")

data("lalonde", package = "MatchIt")

f <- treat ~ age + educ + race + married + nodegree + re74 + re75

lalonde_sw <- seq(0.5, 2, length.out = nrow(lalonde))

#The matched sample must satisfy the requested tolerance on every covariate.
#`matchit2cardinality()` standardizes by the SD of the focal group in the *full*
#sample (not the matched sample), so this must use the same basis or it will report
#violations that are only a difference of convention.
expect_tols_met <- function(m, tols, std.tols = TRUE) {
  X <- m$X[, vapply(m$X, is.numeric, logical(1L)), drop = FALSE]
  keep <- m$weights > 0

  diffs <- vapply(X, function(x) {
    abs(mean(x[keep & m$treat == 1L]) - mean(x[keep & m$treat == 0L]))
  }, numeric(1L))

  if (std.tols) {
    sds <- vapply(X, function(x) sd(x[m$treat == 1L]), numeric(1L))
    diffs <- diffs / ifelse(sds > 0, sds, 1)
  }

  #A little slack for the solver's own tolerance
  expect_lte(max(diffs), tols + 1e-6)

  invisible(m)
}

test_that("baseline: default tols, highs solver", {
  m <- matchit(f, data = lalonde, method = "cardinality")
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = FALSE,
                      expect_subclass = FALSE)
  expect_matchit_snapshot(m)
})

test_that("the requested tolerance is actually met", {
  m <- matchit(f, data = lalonde, method = "cardinality", tols = 0.05)
  expect_tols_met(m, 0.05)
})

test_that("tols = 0.1", {
  m <- matchit(f, data = lalonde, method = "cardinality", tols = 0.1)
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = FALSE,
                      expect_subclass = FALSE)
  expect_tols_met(m, 0.1)
  expect_matchit_snapshot(m)
})

test_that("looser tolerances retain at least as many units", {
  n_kept <- vapply(c(0.01, 0.05, 0.1, 0.25), function(tol) {
    sum(matchit(f, data = lalonde, method = "cardinality", tols = tol)$weights > 0)
  }, numeric(1L))

  expect_false(is.unsorted(n_kept))
})

test_that("std.tols = FALSE", {
  m <- matchit(f, data = lalonde, method = "cardinality",
               tols = 0.1, std.tols = FALSE)
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = FALSE,
                      expect_subclass = FALSE)
  expect_matchit_snapshot(m)
})

test_that("ratio = 2", {
  m <- matchit(f, data = lalonde, method = "cardinality", ratio = 2)
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = FALSE,
                      expect_subclass = FALSE)

  keep <- m$weights > 0
  expect_equal(sum(keep & m$treat == 0L), 2L * sum(keep & m$treat == 1L))

  expect_matchit_snapshot(m)
})

test_that("ratio = NA performs profile matching", {
  m <- matchit(f, data = lalonde, method = "cardinality", ratio = NA)
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = FALSE,
                      expect_subclass = FALSE)

  #Profile matching retains all focal units and selects controls to match them
  expect_true(all(m$weights[m$treat == 1L] > 0))

  expect_matchit_snapshot(m)
})

test_that("estimand='ATC'", {
  m <- matchit(f, data = lalonde, method = "cardinality", estimand = "ATC")
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = FALSE,
                      expect_subclass = FALSE)
  expect_matchit_snapshot(m)
})

test_that("estimand='ATE'", {
  m <- matchit(f, data = lalonde, method = "cardinality", estimand = "ATE")
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = FALSE,
                      expect_subclass = FALSE)
  expect_matchit_snapshot(m)
})

test_that("exact", {
  m <- matchit(f, data = lalonde, method = "cardinality", exact = ~ race)
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = FALSE,
                      expect_subclass = FALSE)
  expect_matchit_snapshot(m)
})

test_that("mahvars pairs units after selection", {
  skip_if_not_installed("optmatch")

  #With `mahvars`, cardinality matching selects units and then pairs them, so
  #unlike every other cardinality specification this one returns a match.matrix
  #and a subclass.
  m <- matchit(f, data = lalonde, method = "cardinality",
               mahvars = ~ age + educ, ratio = 1)
  expect_false(is_null(m$match.matrix))
  expect_false(is_null(m$subclass))
  expect_matchit_snapshot(m)
})

test_that("mahvars with a non-integer ratio is an error", {
  expect_err(matchit(f, data = lalonde, method = "cardinality",
                     mahvars = ~ age, ratio = NA),
             "can only be used with")
})

test_that("solver = 'glpk' agrees with highs", {
  skip_if_not_installed("Rglpk")

  m_highs <- matchit(f, data = lalonde, method = "cardinality", solver = "highs")
  m_glpk <- matchit(f, data = lalonde, method = "cardinality", solver = "glpk")

  #Both solvers should find an optimum of the same size, though not necessarily
  #the same set of units
  expect_equal(sum(m_glpk$weights > 0), sum(m_highs$weights > 0))
  expect_good_matchit(m_glpk, expect_distance = FALSE, expect_match.matrix = FALSE,
                      expect_subclass = FALSE)
})

test_that("solver = 'symphony' is no longer available", {
  #SYMPHONY was removed because it returned a different solution on every call and
  #did not respond to set.seed().
  expect_err(matchit(f, data = lalonde, method = "cardinality",
                     solver = "symphony"),
             '`solver` should be one of "highs", "glpk", or "gurobi".')
})

test_that("matching improves balance", {
  m <- matchit(f, data = lalonde, method = "cardinality")
  expect_balance_improved(m)
})

test_that("tols of the wrong length is an error", {
  expect_err(matchit(f, data = lalonde, method = "cardinality",
                     tols = c(0.1, 0.2)),
             "`tols` must have length equal to 1 or the number of covariates.")
})

test_that("unused arguments warn and are ignored", {
  expect_wrn(
    matchit(f, data = lalonde, method = "cardinality", caliper = 0.1),
    'The argument `caliper` is not used with `method = "cardinality"` and will be ignored.'
  )
})

test_that("non-constant s.weights currently make the problem unsolvable", {
  #KNOWN BUG/PERFORMANCE DEFECT. With constant `s.weights` the problem solves in
  #well under a second; with any non-constant `s.weights` -- even values as mild as
  #runif(0.9, 1.1) -- the solver hits its time limit and errors, at any `tols`
  #tried up to 0.5. `time` is set low here so the test does not spend the default
  #two minutes proving it. See _dev/method-tests-findings.md.
  expect_no_error(matchit(f, data = lalonde, method = "cardinality",
                          s.weights = rep(2, nrow(lalonde))))

  expect_err(matchit(f, data = lalonde, method = "cardinality",
                     s.weights = lalonde_sw, time = 2),
             "failed to find an optimal solution")
})

test_that("no unexpected conditions in the baseline call", {
  expect_no_unexpected_warning(matchit(f, data = lalonde, method = "cardinality"))
})
