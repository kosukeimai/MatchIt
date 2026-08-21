# Tests for print.matchit().
#
# The printed output is the first thing a user sees, and every line of it is assembled
# from `info`, so it is a compact check that `info` records what the call asked for.
# The text is pinned with `expect_snapshot()`, which testthat skips on CRAN, and the
# lines that carry information the object exposes elsewhere are also checked directly
# against that component, so a change in either one is caught without the snapshot.

data("lalonde", package = "MatchIt")

f_pr <- treat ~ age + educ + race + re74

#Capture printed output as a character vector of lines.
print_lines <- function(x, ...) {
  capture.output(print(x, ...))
}

test_that("print.matchit: returns its input invisibly", {
  m <- matchit(f_pr, data = lalonde)

  out <- NULL
  txt <- capture.output(out <- withVisible(print(m)))

  #Printed something, returned the object itself, and returned it invisibly
  expect_gt(length(txt), 1L)
  expect_false(out$visible)
  expect_identical(out$value, m)
})

test_that("print.matchit: the header and the standard lines are present", {
  m <- matchit(f_pr, data = lalonde)
  lines <- print_lines(m)

  expect_identical(lines[1L], "A `matchit` object")

  #Every line after the header is a top-level bullet or a continuation of one
  bullets <- grep("^ - ", lines, value = TRUE)
  expect_true(all(grepl("^( - |    )", lines[-1L])))

  expect_match(bullets, "^ - method: 1:1 nearest neighbor matching without replacement",
               all = FALSE)
  expect_match(bullets, "^ - distance: Propensity score", all = FALSE)
  expect_match(bullets, "^ - target estimand: ATT", all = FALSE)
  expect_match(bullets, "^ - covariates: age, educ, race, re74", all = FALSE)
})

test_that("print.matchit: the sample-size line agrees with the weights", {
  m <- matchit(f_pr, data = lalonde)

  expect_match(print_lines(m), sprintf(" - number of obs.: 614 (original), %s (matched)",
                                       sum(m$weights > 0)),
               all = FALSE, fixed = TRUE)

  #Without matching every unit is retained, so only the original count is printed
  m0 <- matchit(f_pr, data = lalonde, method = NULL)
  expect_match(print_lines(m0), " - number of obs.: 614 (original)", all = FALSE,
               fixed = TRUE)
  expect_no_match(print_lines(m0), "(matched)", fixed = TRUE)
})

test_that("print.matchit: the distance line ends before the next bullet", {
  #Two branches of the distance block used to leave the line unterminated, so the
  #next bullet was printed onto the end of it
  fits <- list(
    mahalanobis = matchit(f_pr, data = lalonde, distance = "robust_mahalanobis"),
    user = matchit(f_pr, data = lalonde, distance = lalonde$re75 / max(lalonde$re75)),
    mahvars = matchit(f_pr, data = lalonde, mahvars = ~ age + educ),
    ps = matchit(f_pr, data = lalonde)
  )

  for (i in names(fits)) {
    lines <- print_lines(fits[[i]])

    #Exactly one line begins each bullet, and no line carries two of them
    expect_length(grep("- number of obs.:", lines, fixed = TRUE), 1L)
    expect_false(any(grepl(" - distance:.* - number of obs.:", lines)), info = i)

    #No blank lines
    expect_false(any(!nzchar(lines)), info = i)
  }
})

test_that("print.matchit: the caliper line names each caliper and its width", {
  m <- matchit(f_pr, data = lalonde, caliper = c(0.1, age = 2))
  lines <- grep("caliper", print_lines(m), value = TRUE)

  #The unnamed distance caliper is labeled and both are shown on the scale used
  expect_match(lines, "<distance>", all = FALSE, fixed = TRUE)
  expect_match(lines, sprintf("age (%s)", format(round(m$caliper["age"], 3L))),
               all = FALSE, fixed = TRUE)

  #The bracketed annotation on the distance line lists every use of the distance
  expect_match(print_lines(m), "Propensity score [matching, caliper]", all = FALSE,
               fixed = TRUE)

  #Matching on the Mahalanobis distance moves "matching" onto its own line, leaving
  #the propensity score with only the uses that remain
  m2 <- matchit(f_pr, data = lalonde, mahvars = ~ age, caliper = 0.1,
                discard = "control")
  expect_match(print_lines(m2), "Mahalanobis [matching]", all = FALSE, fixed = TRUE)
  expect_match(print_lines(m2), "Propensity score [caliper, common support]",
               all = FALSE, fixed = TRUE)

  #Subclassification is named as such, and with no matching neither label appears
  m3 <- matchit(f_pr, data = lalonde, method = "subclass", subclass = 4,
                discard = "both")
  expect_match(print_lines(m3), "Propensity score [subclassification, common support]",
               all = FALSE, fixed = TRUE)

  m4 <- matchit(f_pr, data = lalonde, method = NULL, discard = "control")
  expect_match(print_lines(m4), "Propensity score [common support]", all = FALSE,
               fixed = TRUE)
})

test_that("print.matchit: common support and sampling weights are reported", {
  m <- matchit(f_pr, data = lalonde, discard = "control")
  expect_match(print_lines(m), " - common support: control units dropped", all = FALSE,
               fixed = TRUE)

  m2 <- matchit(f_pr, data = lalonde, discard = "both")
  expect_match(print_lines(m2), "units from both groups", all = FALSE, fixed = TRUE)

  #Sampling weights get their own line, and the distance block records whether they
  #were used in fitting the propensity score
  sw <- rep(c(1, 2), length.out = nrow(lalonde))
  m3 <- matchit(f_pr, data = lalonde, s.weights = sw)
  expect_match(print_lines(m3), " - sampling weights: present", all = FALSE, fixed = TRUE)
  expect_match(print_lines(m3), "sampling weights included in estimation", all = FALSE,
               fixed = TRUE)

  m4 <- add_s.weights(matchit(f_pr, data = lalonde), sw)
  expect_match(print_lines(m4), "sampling weights not included in estimation",
               all = FALSE, fixed = TRUE)

  #No sampling weights, no line
  expect_no_match(print_lines(m), "sampling weights", fixed = TRUE)
})

test_that("print.matchit: the method line describes the method that was used", {
  expected <- c(
    "1:1 nearest neighbor matching without replacement",
    "2:1 nearest neighbor matching with replacement",
    "Exact matching",
    "Coarsened exact matching",
    "Subclassification (4 subclasses)",
    "Cardinality matching",
    "None (no matching)"
  )

  fits <- list(
    matchit(f_pr, data = lalonde),
    matchit(f_pr, data = lalonde, replace = TRUE, ratio = 2),
    matchit(treat ~ race + married, data = lalonde, method = "exact"),
    matchit(f_pr, data = lalonde, method = "cem"),
    matchit(f_pr, data = lalonde, method = "subclass", subclass = 4),
    matchit(f_pr, data = lalonde, method = "cardinality"),
    matchit(f_pr, data = lalonde, method = NULL)
  )

  for (i in seq_along(fits)) {
    expect_match(print_lines(fits[[i]]), sprintf(" - method: %s", expected[i]),
                 all = FALSE, fixed = TRUE)
  }
})

test_that("print.matchit: many covariates are not named", {
  #The cutoff is 40 columns of `X`, counted before dummy expansion
  set.seed(1234)
  d <- cbind(lalonde["treat"],
             as.data.frame(matrix(rnorm(nrow(lalonde) * 41L), nrow(lalonde))))

  m <- matchit(reformulate(setdiff(names(d), "treat"), "treat"), data = d,
               method = NULL)

  expect_match(print_lines(m), " - covariates: too many to name", all = FALSE,
               fixed = TRUE)

  #One fewer and they are all named
  m2 <- matchit(reformulate(setdiff(names(d), c("treat", "V41")), "treat"), data = d,
                method = NULL)
  expect_match(print_lines(m2), " - covariates: V1, V2", all = FALSE)
})

test_that("print.matchit: printed output is stable", {
  skip_on_cran()

  #A propensity score fit exercises the whole distance block; the others each turn on
  #a different part of it
  expect_snapshot(print(matchit(f_pr, data = lalonde)))

  expect_snapshot(print(matchit(f_pr, data = lalonde, method = NULL)))

  expect_snapshot(print(matchit(f_pr, data = lalonde, mahvars = ~ age + educ,
                                caliper = c(0.1, age = 2), discard = "control",
                                replace = TRUE, ratio = 2)))

  expect_snapshot(print(matchit(treat ~ race + married, data = lalonde,
                                method = "exact")))

  expect_snapshot(print(matchit(f_pr, data = lalonde, method = "subclass",
                                subclass = 4)))

  expect_snapshot(print(matchit(f_pr, data = lalonde, method = "cardinality")))

  expect_snapshot(print(matchit(f_pr, data = lalonde,
                                distance = "robust_mahalanobis")))

  expect_snapshot(print(matchit(f_pr, data = lalonde,
                                s.weights = rep(c(1, 2), length.out = nrow(lalonde)))))
})

test_that("print.matchit: printing emits no conditions", {
  fits <- list(
    matchit(f_pr, data = lalonde),
    matchit(f_pr, data = lalonde, method = "cem"),
    matchit(f_pr, data = lalonde, method = "subclass", subclass = 4),
    matchit(f_pr, data = lalonde, method = NULL)
  )

  for (m in fits) {
    expect_silent(invisible(capture.output(print(m))))
  }
})
