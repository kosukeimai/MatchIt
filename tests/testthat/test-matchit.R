test_that("Argument processing works", {
  data("lalonde", package = "MatchIt")

  #Check ignored arguments
  expect_warning({
    m <- matchit(treat ~ age + educ + married, data = lalonde,
                 method = "cem",
                 distance = "probit")
  }, .w("`distance` ignored"))

  #No warning for ignored args set to default
  expect_no_condition({
    m <- matchit(treat ~ age + educ + married, data = lalonde,
                 method = "cem",
                 distance = "glm")
  })

  #Check error arguments
  expect_error({
    m <- matchit(treat ~ age + educ + married, data = lalonde,
                 method = "nearest",
                 distance = "scaled_euclidean",
                 reestimate = TRUE)
  }, .w("`reestimate` not used"))

  #No error for ignored args set to default
  expect_no_condition({
    m <- matchit(treat ~ age + educ + married, data = lalonde,
                 method = "nearest",
                 distance = "scaled_euclidean",
                 reestimate = FALSE)
  })


})
