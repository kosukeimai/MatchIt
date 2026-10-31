# Tests for method = "optimal".
# Structural checks plus snapshot pins recording present behavior. The matching is
# delegated to optmatch, so the pinned values are a joint property of MatchIt and the
# installed optmatch; a snapshot failure after an optmatch upgrade is expected rather
# than a MatchIt regression.

skip_if_not_installed("optmatch")

data("lalonde", package = "MatchIt")

f <- treat ~ age + educ + race + married + nodegree + re74 + re75
fm <- treat ~ age + educ + re74 + re75

lalonde_sw <- seq(0.5, 2, length.out = nrow(lalonde))

test_that("baseline: 1:1 optimal pair matching", {
  m <- matchit(f, data = lalonde, method = "optimal")
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L, replace = FALSE)
  expect_matchit_snapshot(m)
})

test_that("ratio = 2", {
  m <- matchit(f, data = lalonde, method = "optimal", ratio = 2)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 2L, replace = FALSE)
  expect_matchit_snapshot(m)
})

test_that("variable ratio", {
  m <- matchit(f, data = lalonde, method = "optimal",
               ratio = 2, min.controls = 1, max.controls = 4)
  ratio_attr <- structure(2L, min.controls = 1, max.controls = 4)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = ratio_attr, replace = FALSE)
  expect_matchit_snapshot(m)
})

test_that("exact", {
  #`expect_wrn()` lets the call finish, so `m` is available afterward and the match
  #does not have to be run a second time under `suppressWarnings()`
  expect_wrn(
    m <- matchit(f, data = lalonde, method = "optimal", exact = ~ race),
    "Fewer control units than treated units in some `exact` strata; not all treated units will get a match."
  )

  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L, replace = FALSE)
  expect_matchit_snapshot(m)
})

test_that("exact + ratio = 2", {
  m <- suppressWarnings(matchit(f, data = lalonde, method = "optimal",
                               exact = ~ race, ratio = 2))
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 2L, replace = FALSE)
  expect_matchit_snapshot(m)
})

test_that("mahvars", {
  m <- matchit(f, data = lalonde, method = "optimal",
               mahvars = ~ age + educ + re74 + re75)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L, replace = FALSE)
  expect_matchit_snapshot(m)
})

test_that("full Mahalanobis", {
  m <- matchit(fm, data = lalonde, method = "optimal", distance = "mahalanobis")
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L, replace = FALSE)
  expect_matchit_snapshot(m)
})

test_that("estimand='ATC'", {
  expect_wrn(
    m <- matchit(f, data = lalonde, method = "optimal", estimand = "ATC"),
    "Fewer treated units than control units; not all control units will get a match."
  )

  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L, replace = FALSE)
  expect_matchit_snapshot(m)
})

test_that("discard", {
  m <- matchit(f, data = lalonde, method = "optimal", discard = "both")
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L, replace = FALSE)
  expect_matchit_snapshot(m)
})

test_that("exact strata are taken from the distance matrix in one pass", {
  #Each stratum's block built from .infsm_entries_by_stratum() must be identical to
  #subsetting the whole InfinitySparseMatrix, including after a caliper has removed
  #some of its entries
  set.seed(2)
  n1 <- 30
  n0 <- 50
  z <- setNames(c(rep(1, n1), rep(0, n0)), paste0("u", seq_len(n1 + n0)))
  d <- matrix(runif(n1 * n0), n1, n0,
              dimnames = list(names(z)[z == 1], names(z)[z == 0]))

  mo <- optmatch::as.InfinitySparseMatrix(optmatch::match_on(d))
  mo <- mo + optmatch::caliper(mo, .5)

  row_ex <- sample(1:4, n1, TRUE)
  col_ex <- sample(1:4, n0, TRUE)
  entries <- .infsm_entries_by_stratum(mo, row_ex, col_ex, 4L)

  for (e in 1:4) {
    expect_identical(.subset_infsm(mo, entries[[e]])[which(row_ex == e), which(col_ex == e)],
                     mo[row_ex == e, col_ex == e])
  }
})

test_that("discard leaving an exact stratum with one treated and one control unit", {
  #Once u1 and u2 are discarded, stratum "a" holds one treated and one control unit,
  #which are paired directly rather than by optmatch
  d <- data.frame(treat = c(1, 0, 1, 0, 1, 1, 0, 0, 0),
                  x = c(.1, .2, .5, .55, .3, .35, .4, .32, .37),
                  g = c("a", "a", "a", "a", "b", "b", "b", "b", "b"),
                  row.names = paste0("u", 1:9))
  discard <- rownames(d) %in% c("u1", "u2")

  m <- matchit(treat ~ x, data = d, method = "optimal", distance = d$x,
               exact = ~g, discard = discard)

  expect_true(all(is.na(m$subclass[c("u1", "u2")])))
  expect_false(anyNA(m$subclass[c("u3", "u4")]))
  expect_identical(m$subclass[["u3"]], m$subclass[["u4"]])
  expect_identical(unname(m$match.matrix["u3", 1L]), "u4")
})

test_that("s.weights", {
  m <- matchit(f, data = lalonde, method = "optimal", s.weights = lalonde_sw)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = TRUE,
                      expect_subclass = TRUE, ratio = 1L, replace = FALSE)
  expect_matchit_snapshot(m)
})

test_that("include.obj returns the optmatch object", {
  m <- matchit(f, data = lalonde, method = "optimal", include.obj = TRUE)
  expect_s3_class(m$obj, "optmatch")
})

test_that("optimal matching is no worse than greedy on total distance", {
  #The defining property of optimal matching: it minimizes the sum of within-pair
  #distances, so it cannot do worse than greedy matching on that criterion.
  #
  #This only holds once optmatch's `tol` is tightened. At the default (`tol = 1e-3`,
  #which MatchIt does not override) optimal matching returns a solution that is
  #measurably *worse* than greedy on this data -- 39.69374 against 39.69280 -- and
  #at `tol = 1e-6` the two agree to every printed digit. See
  #_dev/method-tests-findings.md.
  #
  #The comparison allows `1e-8` of slack because the two totals can then differ in
  #the last bits: each call fits its own propensity score model, and some BLAS
  #libraries (e.g., MKL) are not bitwise reproducible across calls.
  total_dist <- function(m) {
    sum(vapply(rownames(m$match.matrix), function(t1) {
      t0 <- m$match.matrix[t1, ]
      sum(abs(m$distance[t1] - m$distance[na.omit(t0)]))
    }, numeric(1L)))
  }

  mn <- matchit(f, data = lalonde, method = "nearest")

  mo_tight <- matchit(f, data = lalonde, method = "optimal", tol = 1e-6)
  expect_lte(total_dist(mo_tight), total_dist(mn) + 1e-8)

  mo_default <- matchit(f, data = lalonde, method = "optimal")
  expect_gt(total_dist(mo_default), total_dist(mo_tight))
})

test_that("matching improves balance", {
  m <- matchit(f, data = lalonde, method = "optimal")
  expect_balance_improved(m)
})

test_that("calipers are not supported and warn", {
  expect_wrn(
    m1 <- matchit(f, data = lalonde, method = "optimal", caliper = 0.2),
    'The argument `caliper` is not used with `method = "optimal"` and will be ignored.'
  )

  #Ignoring the caliper must give the same answer as not supplying one
  m0 <- matchit(f, data = lalonde, method = "optimal")
  expect_identical(m0$match.matrix, m1$match.matrix)
})

test_that("mahvars with a full-distance `distance` is an error", {
  expect_err(matchit(fm, data = lalonde, method = "optimal",
                     distance = "mahalanobis",
                     mahvars = ~ age + educ),
             "cannot be used with")
})

test_that("missing values in covariates are an error", {
  lalonde_na <- inject_missingness(lalonde, "educ")

  expect_err(matchit(f, data = lalonde_na, method = "optimal"),
             "Missing and non-finite values are not allowed in the covariates")
})

test_that("antiexact fails informatively when optmatch finds no match", {
  #optmatch reports "Matching failed. (Restrictions impossible to meet?)" for this
  #specification and returns all NA. MatchIt now catches that and errors clearly
  #rather than carrying the empty result into `subclass2mmC()`, which used to fail
  #with "negative length vectors are not allowed".
  #`expect_err()` muffles the warnings optmatch emits on the way, so they need no
  #separate `suppressWarnings()`
  expect_err(
    matchit(f, data = lalonde, method = "optimal", antiexact = ~ married),
    "No matches were found."
  )

  #The same antiexact specification succeeds with method = "full", so the
  #restriction is not inherently unsatisfiable
  m <- matchit(f, data = lalonde, method = "full", antiexact = ~ married)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)
})

test_that("no unexpected conditions in the baseline call", {
  expect_no_unexpected_warning(matchit(f, data = lalonde, method = "optimal"))
})
