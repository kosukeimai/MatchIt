skip()

data("lalonde")
n <- nrow(lalonde)
n1 <- sum(lalonde$treat == 1)
n0 <- n - n1

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

test_that("1:1 NN PSM w/o replacement w/ exact matching works", {
  expect_warning({
    m <- matchit(treat ~ age + educ + race + married + nodegree +
                   re74 + re75, data = lalonde, method = "nearest",
                 ratio = 1, replace = FALSE, distance = "glm",
                 estimand = "ATT",
                 exact = ~race + married)
  }, "Fewer control units than treated units in some `exact` strata; not all treated units will get a match.")

  expect_good_matchit(m, expect_subclass = TRUE, expect_distance = TRUE,
                      expect_match.matrix = TRUE, ratio = 1, replace = FALSE)

  #NAs in MM
  expect_true(anyNA(m$match.matrix))

  expect_false(is.null(m$exact))

  mf_exact <- model.frame(m$exact, data = lalonde)

  for (r in 1:ncol(m$match.matrix)) {
    for (i in 1:nrow(m$match.matrix)) {
      if (!is.na(m$match.matrix[i, r])) {
        for (j in 1:ncol(mf_exact)) {
          expect_equal(mf_exact[rownames(m$match.matrix)[i], j],
                       mf_exact[m$match.matrix[i, r], j])
        }
      }
    }
  }

})

test_that("1:1 NN PSM w/ replacement w/ exact matching works", {
  expect_no_condition({
    m <- matchit(treat ~ age + educ + race + married + nodegree +
                   re74 + re75, data = lalonde, method = "nearest",
                 ratio = 1, replace = TRUE, distance = "glm",
                 estimand = "ATT",
                 exact = ~race + married)
  })

  expect_good_matchit(m, expect_subclass = TRUE, expect_distance = TRUE,
                      expect_match.matrix = TRUE, ratio = 1, replace = TRUE)

  expect_false(is.null(m$exact))

  mf_exact <- model.frame(m$exact, data = lalonde)

  for (r in 1:ncol(m$match.matrix)) {
    for (i in 1:nrow(m$match.matrix)) {
      if (!is.na(m$match.matrix[i, r])) {
        for (j in 1:ncol(mf_exact)) {
          expect_equal(mf_exact[rownames(m$match.matrix)[i], j],
                       mf_exact[m$match.matrix[i, r], j])
        }
      }
    }
  }

})

test_that("1:1 NN PSM w/o replacement w/ anti-exact matching works", {
  expect_no_condition({
    m <- matchit(treat ~ age + educ + race + married + nodegree +
                   re74 + re75, data = lalonde, method = "nearest",
                 ratio = 1, replace = FALSE, distance = "glm",
                 estimand = "ATT",
                 antiexact = ~race + married)
  })

  expect_good_matchit(m, expect_subclass = TRUE, expect_distance = TRUE,
                      expect_match.matrix = TRUE, ratio = 1, replace = FALSE)

  #NAs in MM
  expect_true(anyNA(m$match.matrix))

  expect_false(is.null(m$antiexact))

  mf_antiexact <- model.frame(m$antiexact, data = lalonde)

  #Add test
})

test_that("1:1 NN PSM w/ replacement w/ exact and anti-exact matching works", {
  expect_no_condition({
    m <- matchit(treat ~ age + educ + race + married + nodegree +
                   re74 + re75, data = lalonde, method = "nearest",
                 ratio = 1, distance = "glm",
                 estimand = "ATT", exact = ~race + married, antiexact = ~race + married)
  }, "No units were matched.")


  expect_no_condition({
    m <- matchit(treat ~ age + educ + race + married + nodegree +
                   re74 + re75, data = lalonde, method = "nearest",
                 ratio = 1, distance = "glm",
                 estimand = "ATT", exact = ~race + married, antiexact = ~race + age)
  }, "No units were matched.")


})

test_that("Ratio checks work", {
  expect_error({
    m <- matchit(treat ~ age + educ + race + married + nodegree +
                   re74 + re75, data = lalonde, method = "nearest",
                 ratio = NA, replace = FALSE, distance = "glm",
                 estimand = "ATT")
  }, "`ratio` cannot be `NA`.")
  expect_error({
    m <- matchit(treat ~ age + educ + race + married + nodegree +
                   re74 + re75, data = lalonde, method = "nearest",
                 ratio = NA, replace = FALSE, distance = "glm",
                 estimand = "ATT")
  }, "`ratio` cannot be `NA`.")
  expect_error({
    m <- matchit(treat ~ age + educ + race + married + nodegree +
                   re74 + re75, data = lalonde, method = "nearest",
                 ratio = "hello", replace = FALSE, distance = "glm",
                 estimand = "ATT")
  }, "`ratio` must be a single number greater than or equal to 1.")
  expect_error({
    m <- matchit(treat ~ age + educ + race + married + nodegree +
                   re74 + re75, data = lalonde, method = "nearest",
                 ratio = c(4, 6), replace = FALSE, distance = "glm",
                 estimand = "ATT")
  }, "`ratio` must be a single number greater than or equal to 1.")
  expect_error({
    m <- matchit(treat ~ age + educ + race + married + nodegree +
                   re74 + re75, data = lalonde, method = "nearest",
                 ratio = 0, replace = FALSE, distance = "glm",
                 estimand = "ATT")
  }, "`ratio` must be a single number greater than or equal to 1.")
  expect_error({
    m <- matchit(treat ~ age + educ + race + married + nodegree +
                   re74 + re75, data = lalonde, method = "nearest",
                 ratio = 1.5, replace = FALSE, distance = "glm",
                 estimand = "ATT", max.controls = NULL)
  }, "`ratio` must be a whole number when `max.controls` is not specified.")
  expect_error({
    m <- matchit(treat ~ age + educ + race + married + nodegree +
                   re74 + re75, data = lalonde, method = "nearest",
                 ratio = 1, replace = FALSE, distance = "glm",
                 estimand = "ATT", max.controls = NA)
  }, "`max.controls` must be a single positive number.")
  expect_error({
    m <- matchit(treat ~ age + educ + race + married + nodegree +
                   re74 + re75, data = lalonde, method = "nearest",
                 ratio = 1, replace = FALSE, distance = "glm",
                 estimand = "ATT", max.controls = c("4", "6"))
  }, "`max.controls` must be a single positive number.")
  expect_error({
    m <- matchit(treat ~ age + educ + race + married + nodegree +
                   re74 + re75, data = lalonde, method = "nearest",
                 ratio = 1, replace = FALSE, distance = "glm",
                 estimand = "ATT", max.controls = "hello")
  }, "`max.controls` must be a single positive number.")
  expect_error({
    m <- matchit(treat ~ age + educ + race + married + nodegree +
                   re74 + re75, data = lalonde, method = "nearest",
                 ratio = 1, replace = FALSE, distance = "glm",
                 estimand = "ATT", max.controls = c(1, 5))
  }, "`max.controls` must be a single positive number.")
  expect_error({
    m <- matchit(treat ~ age + educ + race + married + nodegree +
                   re74 + re75, data = lalonde, method = "nearest",
                 ratio = 1, replace = FALSE, distance = "glm",
                 estimand = "ATT", max.controls = 1)
  }, "`ratio` must be greater than 1 for variable ratio matching.")
  expect_error({
    m <- matchit(treat ~ age + educ + race + married + nodegree +
                   re74 + re75, data = lalonde, method = "nearest",
                 ratio = 2, replace = FALSE, distance = "glm",
                 estimand = "ATT", max.controls = 1)
  }, "`max.controls` must be greater than `ratio` for variable ratio matching.")
  expect_error({
    m <- matchit(treat ~ age + educ + race + married + nodegree +
                   re74 + re75, data = lalonde, method = "nearest",
                 ratio = 1, replace = FALSE, distance = "glm",
                 estimand = "ATT", min.controls = NA)
  }, "`min.controls` must be a single positive number.")
  expect_error({
    m <- matchit(treat ~ age + educ + race + married + nodegree +
                   re74 + re75, data = lalonde, method = "nearest",
                 ratio = 1, replace = FALSE, distance = "glm",
                 estimand = "ATT", min.controls = c("4", "6"))
  }, "`min.controls` must be a single positive number.")
  expect_error({
    m <- matchit(treat ~ age + educ + race + married + nodegree +
                   re74 + re75, data = lalonde, method = "nearest",
                 ratio = 1, replace = FALSE, distance = "glm",
                 estimand = "ATT", min.controls = "hello")
  }, "`min.controls` must be a single positive number.")
  expect_error({
    m <- matchit(treat ~ age + educ + race + married + nodegree +
                   re74 + re75, data = lalonde, method = "nearest",
                 ratio = 1, replace = FALSE, distance = "glm",
                 estimand = "ATT", min.controls = c(1, 5))
  }, "`min.controls` must be a single positive number.")
  expect_error({
    m <- matchit(treat ~ age + educ + race + married + nodegree +
                   re74 + re75, data = lalonde, method = "nearest",
                 ratio = 2, replace = FALSE, distance = "glm",
                 estimand = "ATT", min.controls = 0.5, max.controls = 3)
  }, "`min.controls` cannot be less than 1 for variable ratio matching.")
  expect_error({
    m <- matchit(treat ~ age + educ + race + married + nodegree +
                   re74 + re75, data = lalonde, method = "nearest",
                 ratio = 2, replace = FALSE, distance = "glm",
                 estimand = "ATT", min.controls = 2, max.controls = 3)
  }, "`min.controls` must be less than `ratio` for variable ratio matching.")
  expect_error({
    m <- matchit(treat ~ age + educ + race + married + nodegree +
                   re74 + re75, data = lalonde, method = "nearest",
                 ratio = 1, replace = FALSE, distance = "GXM",
                 estimand = "ATT")
  }, "The argument supplied to `distance` is not an allowable value. See `help(\"distance\")` for allowable options.")
  expect_error({
    m <- matchit(treat ~ age + educ + race + married + nodegree +
                   re74 + re75, data = lalonde, method = "nearest",
                 ratio = 1, replace = FALSE, distance = list(1),
                 estimand = "ATT")
  }, "`distance` must be a string with the name of the distance measure to be used or a numeric vector or matrix containing distance measures.")
})
test_that("Method accepts abbreviations", {
  expect_no_error({
    m <- matchit(treat ~ age + educ + race + married + nodegree +
                   re74 + re75, data = lalonde, method = "near",
                 ratio = 1, replace = FALSE, distance = "glm",
                 estimand = "ATT")
  })
})
test_that("Identifying s.weights ", {
  expect_error({
    m <- matchit(treat ~ age + educ + race + married + nodegree +
                   re74 + re75, data = lalonde, method = "nearest",
                 ratio = 1, replace = FALSE, distance = "glm",
                 estimand = "ATT", s.weights = "5")
  }, "The name supplied to `s.weights` must be a variable in `data`.")
  expect_error({
    m <- matchit(lalonde$treat ~ lalonde$age, method = "nearest",
                 ratio = 1, replace = FALSE, distance = "glm",
                 estimand = "ATT", s.weights = "age")
  }, "If `s.weights` is specified a string, a data frame containing the named variable must be supplied to `data`.")
  expect_error({
    m <- matchit(treat ~ age + educ + race + married + nodegree +
                   re74 + re75, data = lalonde, method = "nearest",
                 ratio = 1, replace = FALSE, distance = "glm",
                 estimand = "ATT", s.weights = age)
  }, "The name supplied to `s.weights` must be a variable in `data`.")
})
test_that("Caliper is a numeric or atomic vector", {
  expect_error({
    m <- matchit(treat ~ age + educ + race + married + nodegree +
                   re74 + re75, data = lalonde, method = "nearest",
                 ratio = .5, replace = FALSE, distance = "glm",
                 estimand = "ATT", caliper = c("A"))
  }, "`caliper` must be a numeric vector.")
})