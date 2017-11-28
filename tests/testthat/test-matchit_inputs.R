
library(testthat)
library(MatchIt)

data("lalonde")

set.seed(3254L)
temp_df <- data.frame(
  treat= sample.int(2L, 100, replace=TRUE) - 1,
  x1= rnorm(100),
  x2= rnorm(100),
  x3= rnorm(100),
  z1= c(rep(NA, 50), rnorm(50)),
  z2= c(rep(Inf, 50), rnorm(50)),
  z3= c(rep(NaN, 50), rnorm(50)),
  z4= c(rep(-Inf, 50), rnorm(50))
)

context("check matchit inputs")

test_that("matchit correctly identifies errors", {
  # basic type checking
  expect_error(check_data_missingness(treat ~ abc + educ + black + hispan, data= vector(1:10)))
  expect_error(check_data_missingness(treat ~ abc + educ + black + hispan, data= as.list(1:10)))
  expect_error(check_data_missingness(treat ~ abc + educ + black + hispan, data= NULL))
  
  # check that (non-)existant variables are noticed
  expect_error(check_data_missingness(treat ~ abc + educ + black + hispan, data = lalonde))
  expect_error(check_data_missingness(treat ~ educ + black + hispan + abc, data = lalonde))
  expect_silent(check_data_missingness(treat ~ educ + black + hispan, data = lalonde))
  expect_silent(check_data_missingness(treat ~ x1 + x2 + x3, data= temp_df))
  
  # formula failures
  expect_error(check_data_missingness(1 + 2, data = lalonde))
  expect_error(check_data_missingness(3 ~ 2 + 1, data = lalonde)) # data not exist
  
  # missing value errors
  expect_error(check_data_missingness(treat ~ x1 + x2 + x3 + z1, data= temp_df))
  expect_error(check_data_missingness(treat ~ x1 + x2 + x3 + z2, data= temp_df))
  expect_error(check_data_missingness(treat ~ x1 + x2 + x3 + z3, data= temp_df))
  expect_error(check_data_missingness(treat ~ x1 + x2 + x3 + z4, data= temp_df))
  
  l2 <- lalonde
  l2$educ[sample.int(nrow(l2), 100, replace=FALSE)] <- NA
  
  expect_error(check_data_missingness(treat ~ abc + educ + black + hispan, data = l2))
  
})


test_that("matchit correctly allows missing values", {
  l2 <- lalonde
  l2$educ[sample.int(nrow(l2), 100, replace=FALSE)] <- NA
  
  expect_silent(check_data_missingness(treat ~ x1 + x2 + x3, data= temp_df))
  expect_silent(check_data_missingness(treat ~ black + hispan, data = l2))
  expect_silent(check_data_missingness(treat ~ age + black + hispan, data = l2))
  expect_silent(check_data_missingness(treat ~ age + black + married, data = l2))
})