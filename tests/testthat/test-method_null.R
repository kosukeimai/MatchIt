# Tests for method = NULL, which estimates a distance measure and applies `discard`
# but performs no matching. Everything it returns is deterministic and safe to pin.

data("lalonde", package = "MatchIt")

f <- treat ~ age + educ + race + married + nodegree + re74 + re75

test_that("baseline: no matching, all units retained", {
  m <- matchit(f, data = lalonde, method = NULL)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = FALSE)
  expect_null(m$info$method)
  expect_true(all(m$weights == 1))
  expect_matchit_snapshot(m)
})

test_that("discard drops units and zeroes their weights", {
  m <- matchit(f, data = lalonde, method = NULL, discard = "both")
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = FALSE)

  expect_true(any(m$discarded))
  expect_true(all(m$weights[m$discarded] == 0))
  expect_true(all(m$weights[!m$discarded] == 1))

  expect_matchit_snapshot(m)
})

test_that("the propensity score matches the one from a matching method", {
  #`method = NULL` shares the distance-estimation path with every other method, so
  #the scores must agree exactly.
  m0 <- matchit(f, data = lalonde, method = NULL)
  m1 <- matchit(f, data = lalonde, method = "nearest")

  expect_identical(m0$distance, m1$distance)
})

test_that("s.weights are carried through", {
  sw <- seq(0.5, 2, length.out = nrow(lalonde))
  m <- matchit(f, data = lalonde, method = NULL, s.weights = sw)

  #`in_ps` records whether the weights entered the propensity score model
  expect_equal(m$s.weights, sw, ignore_attr = TRUE)
  expect_true(attr(m$s.weights, "in_ps"))

  expect_matchit_snapshot(m)
})

test_that("distance = 'mahalanobis' is accepted but produces no distance", {
  #A full-sample distance measure has nothing to attach to when no matching
  #happens, so `distance` comes back empty while `info` still records the request.
  m <- matchit(f, data = lalonde, method = NULL, distance = "mahalanobis")

  expect_null(m$distance)
  expect_identical(m$info$distance, "mahalanobis")
  expect_true(all(m$weights == 1))
})

test_that("unused arguments warn and are ignored", {
  expect_wrn(
    matchit(f, data = lalonde, method = NULL, ratio = 2),
    "The argument `ratio` is not used with `method = NULL` and will be ignored."
  )
})

test_that("no unexpected conditions in the baseline call", {
  expect_no_unexpected_warning(matchit(f, data = lalonde, method = NULL))
})
