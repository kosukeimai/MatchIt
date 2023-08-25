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

#  for (r in 1:ncol(m$match.matrix)) {
 #   for (i in 1:nrow(m$match.matrix)) {
#      if (!is.na(m$match.matrix[i, r])) {

 #     }
  #  }
  #}

})

test_that("Mahvars is NULL and reestimate is FALSE", {
  expect_error({
    m <- matchit(treat ~ age + educ + race + married + nodegree +
                   re74 + re75, data = lalonde, method = "nearest",
                 ratio = 1, replace = FALSE, distance = "glm",
                 estimand = "ATT",
                 antiexact = ~race + married, reestimate = TRUE)
  })
  expect_error({
    m <- matchit(treat ~ age + educ + race + married + nodegree +
                   re74 + re75, data = lalonde, method = "nearest",
                 ratio = 1, replace = FALSE, distance = "glm",
                 estimand = "ATT",
                 antiexact = ~race + married, mahvars = 1)
  })

})

test_that("Method is optimal, nearest or full when supplied distance matrix", {
  expect_error({
    m <- matchit(treat ~ age + educ + race + married + nodegree +
                   re74 + re75, data = lalonde, method = "quick",
                 ratio = 1, replace = FALSE, distance = "glm",
                 estimand = "ATT",
                 antiexact = ~race + married)
  })
  expect_error({
    m <- matchit(treat ~ age + educ + race + married + nodegree +
                   re74 + re75, data = lalonde, method = "genetic",
                 ratio = 1, replace = FALSE, distance = "glm",
                 estimand = "ATT",
                 antiexact = ~race + married)
  })
  expect_error({
    m <- matchit(treat ~ age + educ + race + married + nodegree +
                   re74 + re75, data = lalonde, method = "exact",
                 ratio = 1, replace = FALSE, distance = "glm",
                 estimand = "ATT",
                 antiexact = ~race + married)
  })
  expect_error({
    m <- matchit(treat ~ age + educ + race + married + nodegree +
                   re74 + re75, data = lalonde, method = "quick",
                 ratio = 1, replace = FALSE, distance = "cem",
                 estimand = "ATT",
                 antiexact = ~race + married)
  })
  expect_error({
    m <- matchit(treat ~ age + educ + race + married + nodegree +
                   re74 + re75, data = lalonde, method = "subclass",
                 ratio = 1, replace = FALSE, distance = "glm",
                 estimand = "ATT",
                 antiexact = ~race + married)
  })
  expect_error({
    m <- matchit(treat ~ age + educ + race + married + nodegree +
                   re74 + re75, data = lalonde, method = "cardinality",
                 ratio = 1, replace = FALSE, distance = "glm",
                 estimand = "ATT",
                 antiexact = ~race + married)
  })
})

  test_that("Ratio is an integer greater than or equal to 1", {
    expect_error({
      m <- matchit(treat ~ age + educ + race + married + nodegree +
                     re74 + re75, data = lalonde, method = "quick",
                   ratio = 0, replace = FALSE, distance = "glm",
                   estimand = "ATT",
                   antiexact = ~race + married)
    })
    expect_error({
      m <- matchit(treat ~ age + educ + race + married + nodegree +
                     re74 + re75, data = lalonde, method = "quick",
                   ratio = NA, replace = FALSE, distance = "glm",
                   estimand = "ATT",
                   antiexact = ~race + married)
    })
    expect_error({
      m <- matchit(treat ~ age + educ + race + married + nodegree +
                     re74 + re75, data = lalonde, method = "quick",
                   ratio = NULL, replace = FALSE, distance = "glm",
                   estimand = "ATT",
                   antiexact = ~race + married)
    })
  })

  test_that("Ratio must be a whole number when max.controls is NULL", {
    expect_error({
      m <- matchit(treat ~ age + educ + race + married + nodegree +
                     re74 + re75, data = lalonde, method = "quick",
                   ratio = .5, replace = FALSE, distance = "glm",
                   estimand = "ATT",
                   antiexact = ~race + married, max.controls = NULL)
    })
  })

  test_that("Max.controls must be a single numeric value and not contain NAs", {
    expect_error({
      m <- matchit(treat ~ age + educ + race + married + nodegree +
                     re74 + re75, data = lalonde, method = "quick",
                   ratio = .5, replace = FALSE, distance = "glm",
                   estimand = "ATT",
                   antiexact = ~race + married, max.controls = NULL)
    })
    expect_error({
      m <- matchit(treat ~ age + educ + race + married + nodegree +
                     re74 + re75, data = lalonde, method = "quick",
                   ratio = .5, replace = FALSE, distance = "glm",
                   estimand = "ATT",
                   antiexact = ~race + married, max.controls = -1)
    })
    expect_error({
      m <- matchit(treat ~ age + educ + race + married + nodegree +
                     re74 + re75, data = lalonde, method = "quick",
                   ratio = .5, replace = FALSE, distance = "glm",
                   estimand = "ATT",
                   antiexact = ~race + married, max.controls = c(6, 7))
    })


  expect_false(is.null(m$max.controls))
  expect_false(anyNA(m$max.controls))
  })

  test_that("Length of caliper is not 0", {
    expect_error({
      m <- matchit(treat ~ age + educ + race + married + nodegree +
                     re74 + re75, data = lalonde, method = "quick",
                   ratio = .5, replace = FALSE, distance = "glm",
                   estimand = "ATT",
                   antiexact = ~race + married, caliper = 0)
    })
  })
    test_that("Caliper is a numeric or atomic vector", {
      expect_error({
        m <- matchit(treat ~ age + educ + race + married + nodegree +
                       re74 + re75, data = lalonde, method = "quick",
                     ratio = .5, replace = FALSE, distance = "glm",
                     estimand = "ATT",
                     antiexact = ~race + married, caliper = c("A", 2))
      })
    })