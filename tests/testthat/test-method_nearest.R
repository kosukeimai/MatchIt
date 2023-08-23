data("lalonde")
n <- nrow(lalonde)
n1 <- sum(lalonde$treat == 1)
n0 <- n - n1

#put what you've done into a different file, move what you've done from helper
#wait til end of day Friday

test_that("1:1 NN PSM w/o replacement works", {
  expect_no_condition({
    m <- matchit(treat ~ age + educ + race + married + nodegree +
                   re74 + re75, data = lalonde, method = "nearest",
                 ratio = 1, replace = FALSE, distance = "glm",
                 estimand = "ATT")
  })

  expect_good_matchit(m, expect_subclass = TRUE, expect_distance = TRUE,
                      expect_match.matrix = TRUE, ratio = 1, replace = FALSE)

  #No NAs in MM
  expect_false(anyNA(m$match.matrix))

  expect_equal(nlevels(m$subclass), n1)
  expect_equal(sum(!is.na(m$subclass)), 2*n1)

  #More tests...
  expect_length(m$match.matrix, 5:7)
  expect 1 value per row
  expect a string
  expect unique value for each row

})

test_that("2:1 NN PSM w/o replacement works", {
  expect_no_condition({
    m <- matchit(treat ~ age + educ + race + married + nodegree +
                   re74 + re75, data = lalonde, method = "nearest",
                 ratio = 2, replace = FALSE, distance = "glm",
                 estimand = "ATT")
  })

  expect_good_matchit(m, expect_subclass = TRUE, expect_distance = TRUE,
                      expect_match.matrix = TRUE, ratio = 2, replace = FALSE)

  #No NAs in MM
  expect_false(anyNA(m$match.matrix))

  expect_equal(nlevels(m$subclass), n1)
  expect_equal(sum(!is.na(m$subclass)), 3*n1)

  #More tests...
})

test_that("1:1 NN PSM w/o replacement for ATC works", {
  expect_warning({
    m <- matchit(treat ~ age + educ + race + married + nodegree +
                   re74 + re75, data = lalonde, method = "nearest",
                 ratio = 1, replace = FALSE, distance = "glm",
                 estimand = "ATC")
  }, "Fewer treated units than control units; not all control units will get a match")

  expect_good_matchit(m, expect_subclass = TRUE, expect_distance = TRUE,
                      expect_match.matrix = TRUE, ratio = 1, replace = FALSE)

  #No NAs in MM
  expect_equal(sum(is.na(m$match.matrix)), n0 - n1)

  expect_equal(nlevels(m$subclass), n1)
  expect_equal(sum(!is.na(m$subclass)), 2*n1)

  #More tests...
})