# Tests for add_s.weights(), which attaches sampling weights to a fitted `matchit`
# object and, for the stratification methods, recomputes the matching weights from
# them.
#
# For `exact` and `cem` the strata depend only on the covariates, so adding sampling
# weights afterward gives exactly what supplying them to `matchit()` would have. For
# `subclass`, `full`, and `quick` the strata depend on the propensity score, which
# `s.weights` changes, so the two routes legitimately differ -- and that gap is the
# reason the function exists, per its own documentation. The checks below are therefore
# split: exact agreement where the strata match, and an
# implementation-independent balance property everywhere.

data("lalonde", package = "MatchIt")

f3 <- treat ~ age + educ + race
f <- treat ~ age + educ + race + married + nodegree + re74 + re75

lalonde_sw <- seq(0.5, 2, length.out = nrow(lalonde))

strat_methods <- c("exact", "cem", "subclass", "full", "quick")

fit <- function(method, ...) {
  matchit(if (method == "exact") f3 else f, data = lalonde, method = method, ...)
}

#Within every stratum the ratio of weighted control mass to weighted treated mass must
#be the same. `matchit()` rescales the nonzero weights to mean 1 within each treatment
#group, which multiplies that ratio by one global constant but cannot make it vary
#across strata, so constancy is the property to check rather than any particular value.
expect_stratum_mass_balanced <- function(m, s.weights) {
  keep <- !is.na(m$subclass) & m$weights > 0

  ratios <- vapply(levels(droplevels(m$subclass[keep])), function(lev) {
    i <- which(keep & m$subclass == lev)
    sum((m$weights * s.weights)[i[m$treat[i] == 0L]]) /
      sum((m$weights * s.weights)[i[m$treat[i] == 1L]])
  }, numeric(1L))

  expect_equal(max(ratios), min(ratios))

  invisible(m)
}

test_that("s.weights are stored with the right in_ps flag", {
  m <- add_s.weights(fit("exact"), lalonde_sw)
  expect_equal(m$s.weights, lalonde_sw, ignore_attr = TRUE)
  expect_false(attr(m$s.weights, "in_ps"))

  #Adding the same weights that were already used in the propensity score
  m <- add_s.weights(fit("exact", s.weights = lalonde_sw), lalonde_sw)
  expect_true(attr(m$s.weights, "in_ps"))
})

test_that("s.weights = NULL returns the object unchanged", {
  m <- fit("exact")
  expect_identical(add_s.weights(m, NULL), m)
})

test_that("the matching weights are recomputed for every stratification method", {
  for (method in strat_methods) {
    m0 <- fit(method)
    m1 <- add_s.weights(m0, lalonde_sw)

    expect_not_equal(m1$weights, m0$weights)

    #The strata themselves are not re-derived
    expect_identical(m1$subclass, m0$subclass)
  }
})

test_that("recomputed weights balance the weighted stratum masses", {
  for (method in strat_methods) {
    expect_stratum_mass_balanced(add_s.weights(fit(method), lalonde_sw),
                                lalonde_sw)
  }
})

test_that("exact and cem agree with supplying s.weights to matchit()", {
  #Their strata do not depend on the propensity score, so the two routes must agree
  for (method in c("exact", "cem")) {
    expect_identical(fit(method)$subclass,
                     fit(method, s.weights = lalonde_sw)$subclass)

    expect_equal(add_s.weights(fit(method), lalonde_sw)$weights,
                 fit(method, s.weights = lalonde_sw)$weights)
  }
})

test_that("PS-based methods differ only because the strata differ", {
  for (method in c("subclass", "full", "quick")) {
    direct <- fit(method, s.weights = lalonde_sw)

    #s.weights changed the propensity score, so the strata are not the same
    expect_not_equal(direct$subclass, fit(method)$subclass)

    #Given the same strata, the two routes agree exactly
    stripped <- direct
    stripped$s.weights <- NULL

    expect_equal(add_s.weights(stripped, lalonde_sw)$weights, direct$weights)
  }
})

test_that("normalize = FALSE from the original call is honored", {
  for (method in strat_methods) {
    direct <- fit(method, s.weights = lalonde_sw, normalize = FALSE)
    expect_false(direct$info$normalize)

    stripped <- direct
    stripped$s.weights <- NULL

    expect_equal(add_s.weights(stripped, lalonde_sw)$weights, direct$weights)
  }
})

test_that("normalize = TRUE and FALSE give different weights", {
  m_norm <- fit("exact", s.weights = lalonde_sw)
  m_raw <- fit("exact", s.weights = lalonde_sw, normalize = FALSE)

  expect_true(m_norm$info$normalize)
  expect_not_equal(m_norm$weights, m_raw$weights)

  #...but the same balance property, since normalization is a per-group rescaling
  expect_stratum_mass_balanced(m_norm, lalonde_sw)
  expect_stratum_mass_balanced(m_raw, lalonde_sw)
})

test_that("normalize is read from info, not from the recorded call", {
  #A call like `matchit(..., normalize = nrm)` records the symbol `nrm`, which cannot
  #be resolved later; `info$normalize` holds the evaluated value instead.
  nrm <- FALSE
  m <- fit("exact", normalize = nrm)

  expect_false(m$info$normalize)
  expect_equal(add_s.weights(m, lalonde_sw)$weights,
               fit("exact", s.weights = lalonde_sw, normalize = FALSE)$weights)
})

test_that("objects predating info$normalize default to normalizing", {
  m <- fit("exact")
  m$info$normalize <- NULL

  expect_equal(add_s.weights(m, lalonde_sw)$weights,
               fit("exact", s.weights = lalonde_sw)$weights)
})

test_that("methods that do not stratify are left alone", {
  for (method in c("nearest", "optimal", "cardinality")) {
    m0 <- suppressWarnings(matchit(f, data = lalonde, method = method))
    m1 <- add_s.weights(m0, lalonde_sw)

    expect_identical(m1$weights, m0$weights)
    expect_equal(m1$s.weights, lalonde_sw, ignore_attr = TRUE)
  }

  m0 <- matchit(f, data = lalonde, method = NULL)
  expect_identical(add_s.weights(m0, lalonde_sw)$weights, m0$weights)
})

test_that("add_s.weights() is idempotent", {
  m1 <- add_s.weights(fit("exact"), lalonde_sw)
  expect_equal(add_s.weights(m1, lalonde_sw)$weights, m1$weights)
})

test_that("nn is recomputed", {
  m0 <- fit("exact")
  m1 <- add_s.weights(m0, lalonde_sw)

  expect_not_equal(m1$nn, m0$nn)
})

test_that("print and summary work on the result", {
  m <- add_s.weights(fit("exact"), lalonde_sw)

  expect_no_error(invisible(capture.output(print(m))))
  expect_no_error(invisible(capture.output(summary(m))))
})

test_that("s.weights can be given as a formula or a name", {
  lalonde_w <- lalonde
  lalonde_w$sw <- lalonde_sw

  target <- add_s.weights(fit("exact"), lalonde_sw)$weights

  expect_equal(add_s.weights(fit("exact"), ~ sw, data = lalonde_w)$weights, target)
  expect_equal(add_s.weights(fit("exact"), "sw", data = lalonde_w)$weights, target)
})

test_that("mismatched length is an error", {
  expect_err(add_s.weights(fit("exact"), lalonde_sw[-1L]),
             "`s.weights` must be the same length as the treatment vector.")
})
