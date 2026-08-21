test_that("Argument processing works", {
  data("lalonde", package = "MatchIt")

  #Check ignored arguments
  expect_wrn({
    m <- matchit(treat ~ age + educ + married, data = lalonde,
                 method = "cem",
                 distance = "probit")
  }, 'The argument `distance` is not used with `method = "cem"` and will be ignored.')

  #No warning for ignored args set to default
  expect_no_condition({
    m <- matchit(treat ~ age + educ + married, data = lalonde,
                 method = "cem",
                 distance = "glm")
  })

  #Check error arguments
  expect_err({
    m <- matchit(treat ~ age + educ + married, data = lalonde,
                 method = "nearest",
                 distance = "scaled_euclidean",
                 reestimate = TRUE)
  }, 'The argument `reestimate` cannot be used with `method = "nearest"` and `distance = "scaled_euclidean"`.')

  #No error for ignored args set to default
  expect_no_condition({
    m <- matchit(treat ~ age + educ + married, data = lalonde,
                 method = "nearest",
                 distance = "scaled_euclidean",
                 reestimate = FALSE)
  })


})
