# Tests for method = "subclass".
# Structural checks plus snapshot pins recording present behavior. Subclassification
# is deterministic given the propensity score, so its results are safe to pin.

data("lalonde", package = "MatchIt")

f <- treat ~ age + educ + race + married + nodegree + re74 + re75

lalonde_sw <- seq(0.5, 2, length.out = nrow(lalonde))

test_that("baseline: default 6 subclasses", {
  m <- matchit(f, data = lalonde, method = "subclass")
  expect_s3_class(m, "matchit.subclass")
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)
  expect_length(levels(m$subclass), 6L)
  expect_matchit_snapshot(m)
  expect_snapshot_value(unname(round(m$q.cut, 8L)), style = "json2")
})

test_that("subclass = 10", {
  m <- matchit(f, data = lalonde, method = "subclass", subclass = 10)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)
  expect_length(levels(m$subclass), 10L)
  expect_matchit_snapshot(m)
  expect_snapshot_value(unname(round(m$q.cut, 8L)), style = "json2")
})

test_that("subclass supplied as a vector of quantiles", {
  m <- matchit(f, data = lalonde, method = "subclass",
               subclass = c(.25, .5, .75))
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)
  expect_length(levels(m$subclass), 4L)
  expect_matchit_snapshot(m)
  expect_snapshot_value(unname(round(m$q.cut, 8L)), style = "json2")
})

test_that("quantiles are padded with 0 and 1 when not supplied", {
  m1 <- matchit(f, data = lalonde, method = "subclass", subclass = c(.25, .5, .75))
  m2 <- matchit(f, data = lalonde, method = "subclass", subclass = c(0, .25, .5, .75, 1))

  expect_identical(m1$subclass, m2$subclass)
  expect_equal(m1$q.cut, m2$q.cut)
})

test_that("min.n = 0 drops strata missing a treatment group", {
  m <- matchit(f, data = lalonde, method = "subclass", subclass = 10, min.n = 0)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)
  expect_matchit_snapshot(m)
})

test_that("min.n = 3 scoots units to fill strata", {
  m <- matchit(f, data = lalonde, method = "subclass", subclass = 10, min.n = 3)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)

  keep <- !is.na(m$subclass)
  tab <- table(m$subclass[keep], m$treat[keep])
  expect_true(all(tab >= 3L))

  expect_matchit_snapshot(m)
})

test_that("min.n = 1 (default) leaves at least one unit per group per stratum", {
  m <- matchit(f, data = lalonde, method = "subclass", subclass = 20)

  keep <- !is.na(m$subclass)
  tab <- table(m$subclass[keep], m$treat[keep])
  expect_true(all(tab >= 1L))
})

test_that("estimand='ATC'", {
  m <- matchit(f, data = lalonde, method = "subclass", estimand = "ATC")
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)
  expect_matchit_snapshot(m)
  expect_snapshot_value(unname(round(m$q.cut, 8L)), style = "json2")
})

test_that("estimand='ATE'", {
  m <- matchit(f, data = lalonde, method = "subclass", estimand = "ATE")
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)
  expect_matchit_snapshot(m)
  expect_snapshot_value(unname(round(m$q.cut, 8L)), style = "json2")
})

test_that("the estimand changes the cut points and the weights", {
  ms <- lapply(c("ATT", "ATC", "ATE"), function(e) {
    matchit(f, data = lalonde, method = "subclass", estimand = e)
  })

  expect_not_equal(ms[[1L]]$q.cut, ms[[2L]]$q.cut)
  expect_not_equal(ms[[1L]]$weights, ms[[2L]]$weights)
  expect_not_equal(ms[[1L]]$weights, ms[[3L]]$weights)
})

test_that("s.weights", {
  m <- matchit(f, data = lalonde, method = "subclass", s.weights = lalonde_sw)
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)
  expect_matchit_snapshot(m)
})

test_that("discard", {
  m <- matchit(f, data = lalonde, method = "subclass", discard = "both")
  expect_good_matchit(m, expect_distance = TRUE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)
  expect_true(any(is.na(m$subclass)))
  expect_matchit_snapshot(m)
})

test_that("all units are retained without discard", {
  m <- matchit(f, data = lalonde, method = "subclass")
  expect_false(anyNA(m$subclass))
  expect_true(all(m$weights > 0))
})

test_that("requesting more subclasses than the PS supports warns", {
  expect_wrn(
    m <- matchit(f, data = lalonde, method = "subclass", subclass = 100),
    "Due to discreteness in the distance measure, fewer subclasses were generated than were requested."
  )

  expect_lt(nlevels(m$subclass), 100L)
  expect_matchit_snapshot(m)
})

test_that("matching improves balance", {
  m <- matchit(f, data = lalonde, method = "subclass", subclass = 6)
  expect_balance_improved(m)
})

test_that("subclass must ask for more than one subclass", {
  #Values that would produce a single stratum containing every unit are rejected
  #rather than silently performing no subclassification.
  bad <- list(1, 0, c(0, 1), 2.5, -1, NA, "a")

  for (s in bad) {
    expect_err(matchit(f, data = lalonde, method = "subclass", subclass = s),
               "`subclass` must either be the number of desired subclasses")
  }
})

test_that("a single quantile is read as a quantile, not a count", {
  #This is the input that used to segfault: it was read as a count, `round(0.5)`
  #gave 0, and the resulting single cut point reached `has_n_unique()` with a
  #non-integer n.
  m <- matchit(f, data = lalonde, method = "subclass", subclass = 0.5)
  expect_length(levels(m$subclass), 2L)

  #Identical to spelling the same quantile vector out in full
  m_full <- matchit(f, data = lalonde, method = "subclass",
                    subclass = c(0, 0.5, 1))
  expect_identical(m$subclass, m_full$subclass)
  expect_equal(m$q.cut, m_full$q.cut)

  expect_matchit_snapshot(m)
})

test_that("duplicated quantiles are collapsed", {
  m1 <- matchit(f, data = lalonde, method = "subclass", subclass = c(.5, .5))
  m2 <- matchit(f, data = lalonde, method = "subclass", subclass = 0.5)

  expect_identical(m1$subclass, m2$subclass)
})

test_that("distance = 'mahalanobis' is an error", {
  expect_err(matchit(f, data = lalonde, method = "subclass",
                     distance = "mahalanobis"),
             '`distance` cannot be "mahalanobis" with `method = "subclass"`.')
})

test_that("sub.by is defunct", {
  expect_err(matchit(f, data = lalonde, method = "subclass", sub.by = "treat"),
             "`sub.by` is defunct and has been replaced with `estimand`.")
})

test_that("unused arguments warn and are ignored", {
  expect_wrn(
    matchit(f, data = lalonde, method = "subclass", caliper = 0.1),
    'The argument `caliper` is not used with `method = "subclass"` and will be ignored.'
  )

  expect_wrn(
    matchit(f, data = lalonde, method = "subclass", exact = ~ race),
    'The argument `exact` is not used with `method = "subclass"` and will be ignored.'
  )
})

test_that("no unexpected conditions in the baseline call", {
  expect_no_unexpected_warning(matchit(f, data = lalonde, method = "subclass"))
})
