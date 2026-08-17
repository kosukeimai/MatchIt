# Tests for method = "exact".
# Structural checks plus snapshot pins recording present behavior. Exact matching is
# fully deterministic and has no external solver, so its results are safe to pin.

data("lalonde", package = "MatchIt")

# Three covariates give 42 populated strata; all seven give only 7, and the full
# formula is used below to exercise the sparse case.
f3 <- treat ~ age + educ + race
f7 <- treat ~ age + educ + race + married + nodegree + re74 + re75

lalonde_sw <- seq(0.5, 2, length.out = nrow(lalonde))

test_that("baseline: three covariates", {
  m <- matchit(f3, data = lalonde, method = "exact")
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)
  expect_matchit_snapshot(m)
})

test_that("all seven covariates (sparse strata)", {
  m <- matchit(f7, data = lalonde, method = "exact")
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)
  expect_matchit_snapshot(m)
})

test_that("single covariate", {
  m <- matchit(treat ~ race, data = lalonde, method = "exact")
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)
  expect_matchit_snapshot(m)
})

test_that("estimand='ATC'", {
  m <- matchit(f3, data = lalonde, method = "exact", estimand = "ATC")
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)
  expect_matchit_snapshot(m)
})

test_that("estimand='ATE'", {
  m <- matchit(f3, data = lalonde, method = "exact", estimand = "ATE")
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)
  expect_matchit_snapshot(m)
})

test_that("s.weights", {
  m <- matchit(f3, data = lalonde, method = "exact", s.weights = lalonde_sw)
  expect_good_matchit(m, expect_distance = FALSE, expect_match.matrix = FALSE,
                      expect_subclass = TRUE)
  expect_matchit_snapshot(m)
})

test_that("the estimand changes the weights but not the strata", {
  ms <- lapply(c("ATT", "ATC", "ATE"), function(e) {
    matchit(f3, data = lalonde, method = "exact", estimand = e)
  })

  expect_identical(ms[[1L]]$subclass, ms[[2L]]$subclass)
  expect_identical(ms[[1L]]$subclass, ms[[3L]]$subclass)

  expect_not_equal(ms[[1L]]$weights, ms[[2L]]$weights)
  expect_not_equal(ms[[1L]]$weights, ms[[3L]]$weights)
})

test_that("s.weights enter the matching weights but not the strata", {
  m0 <- matchit(f3, data = lalonde, method = "exact")
  m1 <- matchit(f3, data = lalonde, method = "exact", s.weights = lalonde_sw)

  expect_identical(m0$subclass, m1$subclass)
  expect_not_equal(m0$weights, m1$weights)
})

test_that("constant s.weights reproduce the unweighted matching weights", {
  m0 <- matchit(f3, data = lalonde, method = "exact")

  for (v in c(1, 3)) {
    m <- matchit(f3, data = lalonde, method = "exact",
                 s.weights = rep(v, nrow(lalonde)))
    expect_equal(m$weights, m0$weights)
  }
})

#Within every stratum the ratio of weighted control mass to weighted treated mass
#must be the same. `matchit()` rescales the nonzero weights to mean 1 within each
#treatment group, which multiplies that ratio by one global constant but cannot make
#it vary across strata -- so constancy is the property to check, not the value.
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

test_that("weighted stratum masses balance, with and without s.weights", {
  expect_stratum_mass_balanced(matchit(f3, data = lalonde, method = "exact"),
                              rep(1, nrow(lalonde)))

  expect_stratum_mass_balanced(
    matchit(f3, data = lalonde, method = "exact", s.weights = lalonde_sw),
    lalonde_sw
  )

  expect_stratum_mass_balanced(
    matchit(f3, data = lalonde, method = "exact", s.weights = lalonde_sw,
            estimand = "ATC"),
    lalonde_sw
  )
})

test_that("add_s.weights() reproduces matching with s.weights", {
  #Exact matching's strata depend only on the covariates, so adding sampling weights
  #after the fact must give exactly what supplying them up front would have.
  #`test-add_s.weights.R` covers the other methods and the `normalize` interaction.
  m0 <- matchit(f3, data = lalonde, method = "exact")
  m1 <- add_s.weights(m0, lalonde_sw)
  m2 <- matchit(f3, data = lalonde, method = "exact", s.weights = lalonde_sw)

  expect_equal(m1$weights, m2$weights)
})

test_that("every retained stratum contains both treatment groups", {
  m <- matchit(f3, data = lalonde, method = "exact")

  keep <- !is.na(m$subclass)
  tab <- table(m$subclass[keep], m$treat[keep])

  expect_true(all(tab[, "0"] > 0L))
  expect_true(all(tab[, "1"] > 0L))
})

test_that("matching improves balance", {
  m <- matchit(f3, data = lalonde, method = "exact")
  expect_balance_improved(m)
})

test_that("no covariates is an error", {
  expect_err(matchit(treat ~ 1, data = lalonde, method = "exact"),
             "Covariates must be specified in the input formula to use exact matching.")
})

test_that("no shared strata is an error", {
  #A covariate unique to every row leaves no value shared across groups. `re74`
  #does not work for this: two thirds of its values are 0, so strata do overlap.
  lalonde_uniq <- lalonde
  lalonde_uniq$uniq <- seq_len(nrow(lalonde))

  expect_err(matchit(treat ~ uniq, data = lalonde_uniq, method = "exact"),
             "No exact matches were found.")
})

test_that("unused arguments warn and are ignored", {
  expect_wrn(
    matchit(f3, data = lalonde, method = "exact", caliper = 0.1),
    'The argument `caliper` is not used with `method = "exact"` and will be ignored.'
  )

  expect_wrn(
    matchit(f3, data = lalonde, method = "exact", replace = TRUE),
    'The argument `replace` is not used with `method = "exact"` and will be ignored.'
  )

  #Supplying both produces one pluralized warning rather than two
  expect_wrn(
    m1 <- matchit(f3, data = lalonde, method = "exact",
                  caliper = 0.1, replace = TRUE),
    'The arguments `caliper` and `replace` are not used with `method = "exact"` and will be ignored.'
  )

  m0 <- matchit(f3, data = lalonde, method = "exact")
  expect_identical(m0$subclass, m1$subclass)
  expect_identical(m0$weights, m1$weights)
})

test_that("missing values in covariates are an error", {
  lalonde_na <- inject_missingness(lalonde, "educ")

  expect_err(matchit(f3, data = lalonde_na, method = "exact"),
             "Missing and non-finite values are not allowed in the covariates")
})

test_that("no unexpected conditions in the baseline call", {
  expect_no_unexpected_warning(matchit(f3, data = lalonde, method = "exact"))
})
