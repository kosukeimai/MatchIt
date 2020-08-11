
library(testthat)
library(MatchIt)
data("lalonde")

context("hist.pscore")

test_that("hist.pscore correctly restores par() settings", {
  # setup test
  m.out <- matchit(treat ~ educ + race, data = lalonde,
                   method = "nearest")

  par_mfrow_prior <- par(no.readonly = TRUE)$mfrow
  plot(m.out, type="hist")
  par_mfrow_now <- par(no.readonly = TRUE)$mfrow

  # tests
  expect_equal(par_mfrow_prior, par_mfrow_now)
})
