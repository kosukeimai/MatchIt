
library(testthat)
library(MatchIt)
data("lalonde")

context("hist.pscore")

test_that("hist.pscore correctly restores par() settings", {
  # setup test
  m.out <- matchit(treat ~ educ + black + hispan, data = lalonde,
                   method = "nearest")
  
  par_prior <- par(no.readonly = TRUE)$mfrow
  plot(m.out, type="hist")
  par_now <- par(no.readonly = TRUE)$mfrow
  
  # tests
  expect_equal(par_prior, par_now)
})
