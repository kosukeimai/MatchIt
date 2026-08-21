# Tests for plot.matchit(), plot.matchit.subclass(), and plot.summary.matchit().
#
# Base graphics leave nothing to inspect, so these tests assert what can be asserted
# without a rendering backend: that each plot type and argument combination runs to
# completion on a null device, that the arguments that select variables and subclasses
# are honored or rejected as documented, and that nothing prompts for input. That last
# one is the point of several of these -- a plot that blocks waiting on `readline()`
# hangs `R CMD check` and knitr rather than failing, so it is worth a test.

data("lalonde", package = "MatchIt")

f_pl <- treat ~ age + educ + race + re74

test_that("plot.matchit: returns its input invisibly", {
  local_null_device()

  m <- matchit(f_pl, data = lalonde)

  out <- withVisible(plot(m, interactive = FALSE))
  expect_false(out$visible)
  expect_identical(out$value, m)

  ms <- matchit(f_pl, data = lalonde, method = "subclass", subclass = 3)
  out_s <- withVisible(plot(ms, interactive = FALSE))
  expect_false(out_s$visible)
  expect_identical(out_s$value, ms)
})

test_that("plot.matchit: every type runs and emits nothing", {
  local_null_device()

  m <- matchit(f_pl, data = lalonde)

  for (type in c("qq", "ecdf", "density", "jitter", "histogram")) {
    expect_silent(plot(m, type = type, interactive = FALSE))
  }

  #`type` is matched partially, as the other arguments of the package are
  expect_silent(plot(m, type = "hist", interactive = FALSE))

  expect_err(plot(m, type = "boxplot", interactive = FALSE), "type")
})

test_that("plot.matchit: jitter and histogram need a distance measure", {
  local_null_device()

  #Mahalanobis matching produces no distance component
  m <- matchit(f_pl, data = lalonde, distance = "mahalanobis")
  expect_null(m$distance)

  for (type in c("jitter", "histogram")) {
    expect_err(plot(m, type = type, interactive = FALSE),
               "cannot be used if no distance measure was estimated or supplied")
  }

  #The covariate plots do not
  expect_silent(plot(m, type = "qq", interactive = FALSE))

  #Same for a method that never estimates one
  mc <- matchit(f_pl, data = lalonde, method = "cem")
  expect_err(plot(mc, type = "jitter", interactive = FALSE), "cannot be used")
})

test_that("plot.matchit: which.xs selects the variables plotted", {
  local_null_device()

  m <- matchit(f_pl, data = lalonde)

  #A character vector and a formula are both accepted
  expect_silent(plot(m, which.xs = c("age", "educ"), interactive = FALSE))
  expect_silent(plot(m, which.xs = ~ age + educ, interactive = FALSE))

  #A formula can transform the variables
  expect_silent(plot(m, which.xs = ~ I(age^2), interactive = FALSE))

  #A name that is nowhere to be found is an error
  expect_err(plot(m, which.xs = "nonesuch", interactive = FALSE),
             "all variables in `which.xs` must be in the supplied <matchit> object or in `data`")

  #Neither a character vector nor a formula
  expect_err(plot(m, which.xs = 1:2, interactive = FALSE),
             "`which.xs` must be supplied as a character vector of names or a one-sided formula")
})

test_that("plot.matchit: data supplies variables the object does not hold", {
  local_null_device()

  m <- matchit(f_pl, data = lalonde)

  lalonde_extra <- lalonde
  lalonde_extra$extra <- lalonde$age * 2

  expect_silent(plot(m, which.xs = "extra", data = lalonde_extra, interactive = FALSE))
  expect_silent(plot(m, which.xs = ~ extra + age, data = lalonde_extra,
                     interactive = FALSE))

  #A wrongly sized `data` is an error rather than being quietly replaced, and so is
  #one that is not a data frame
  expect_err(plot(m, which.xs = "extra", data = lalonde_extra[1:10, ],
                  interactive = FALSE),
             "`data` must have one row for each of the 614 units")
  expect_err(plot(m, which.xs = "extra", data = as.list(lalonde_extra),
                  interactive = FALSE),
             "`data` must be a data frame")

  #Validated even when `which.xs` is not supplied, so a mistake is not silently ignored
  expect_err(plot(m, data = lalonde_extra[1:10, ], interactive = FALSE),
             "`data` must have one row for each of the 614 units")
})

test_that("plot.matchit: works after every matching method", {
  local_null_device()

  fits <- list(
    nearest = matchit(f_pl, data = lalonde),
    replace = matchit(f_pl, data = lalonde, replace = TRUE, ratio = 2),
    cem = matchit(f_pl, data = lalonde, method = "cem"),
    exact = matchit(treat ~ race + married, data = lalonde, method = "exact"),
    cardinality = matchit(f_pl, data = lalonde, method = "cardinality"),
    `no matching` = matchit(f_pl, data = lalonde, method = NULL)
  )

  for (i in names(fits)) {
    for (type in c("qq", "ecdf", "density")) {
      expect_silent(plot(fits[[i]], type = type, interactive = FALSE))
    }
  }

  #Sampling weights are used in the plotted distributions, so the weighted paths run too
  msw <- matchit(f_pl, data = lalonde,
                 s.weights = seq(0.5, 2, length.out = nrow(lalonde)))
  for (type in c("qq", "ecdf", "density")) {
    expect_silent(plot(msw, type = type, interactive = FALSE))
  }
})

test_that("plot.matchit: factor and binary covariates take the categorical path", {
  local_null_device()

  #Factors are plotted as bars rather than points, a separate branch for each type
  m <- matchit(treat ~ race + married + age, data = lalonde)

  for (type in c("qq", "ecdf", "density")) {
    expect_silent(plot(m, type = type, which.xs = ~ race, interactive = FALSE))
    expect_silent(plot(m, type = type, which.xs = c("married", "age"),
                       interactive = FALSE))
  }
})

test_that("plot.matchit: variables used only in exact or mahvars are still available", {
  local_null_device()

  #These are not in `X` for every method, so they are added back before plotting
  expect_wrn(m <- matchit(treat ~ age + educ, data = lalonde, exact = ~ race,
                          mahvars = ~ age + educ),
             "fewer control units than treated units in some `exact` strata")

  expect_silent(plot(m, which.xs = "race", interactive = FALSE))
  expect_silent(plot(m, interactive = FALSE))
})

test_that("plot.matchit.subclass: subclass selects which subclasses are shown", {
  local_null_device()

  ms <- matchit(f_pl, data = lalonde, method = "subclass", subclass = 3)

  expect_silent(plot(ms, subclass = TRUE, interactive = FALSE))
  expect_silent(plot(ms, subclass = FALSE, interactive = FALSE))
  expect_silent(plot(ms, subclass = 1:2, interactive = FALSE))
  expect_silent(plot(ms, subclass = 2, interactive = FALSE))

  #Out of range, and not an index at all
  expect_err(plot(ms, subclass = 99, interactive = FALSE),
             "`subclass` should be TRUE, FALSE, or a vector of subclass indices")
  #A string passes the index check by coercion and is caught one level down instead
  expect_err(plot(ms, subclass = "2", interactive = FALSE),
             "the argument supplied to `subclass` is not the index of any subclass")

  #`which.xs` and `data` work the same way per subclass
  expect_silent(plot(ms, subclass = 1, which.xs = ~ age, interactive = FALSE))
  expect_err(plot(ms, subclass = 1, which.xs = ~ age, data = lalonde[1:10, ],
                  interactive = FALSE),
             "`data` must have one row for each of the 614 units")
})

test_that("plot.matchit.subclass: the interactive menu is not entered non-interactively", {
  #The menu reads from stdin in a loop that only exits on a valid choice, so entering
  #it in a non-interactive session spins forever rather than erroring. This test would
  #hang, not fail, if the guard were removed.
  local_null_device()

  ms <- matchit(f_pl, data = lalonde, method = "subclass", subclass = 3)

  expect_false(interactive())

  #`interactive = TRUE` is the default, and `subclass` is missing, which is the
  #combination that used to reach the menu
  expect_silent(plot(ms))
  expect_silent(plot(ms, type = "ecdf"))

  #The distance-based plots have no menu and are unaffected
  expect_silent(plot(ms, type = "histogram"))
})

test_that("plot.matchit.subclass: jitter and histogram bypass the subclass logic", {
  local_null_device()

  ms <- matchit(f_pl, data = lalonde, method = "subclass", subclass = 3)

  expect_silent(plot(ms, type = "jitter", interactive = FALSE))
  expect_silent(plot(ms, type = "histogram"))

  #A user-supplied distance takes the same path
  expect_wrn(msm <- matchit(f_pl, data = lalonde, method = "subclass", subclass = 3,
                            distance = lalonde$re75 / max(lalonde$re75)),
             "fewer subclasses were generated than were requested")
  expect_silent(plot(msm, type = "jitter", interactive = FALSE))
})

# ===== plot.summary.matchit (the Love plot) =====

test_that("plot.summary.matchit: returns its input invisibly and restores par()", {
  local_null_device()

  s <- summary(matchit(f_pl, data = lalonde))

  before <- par(no.readonly = TRUE)

  out <- withVisible(plot(s))
  expect_false(out$visible)
  expect_identical(out$value, s)

  #The graphical parameters the function changes are put back
  expect_identical(par(no.readonly = TRUE)[names(before)], before)
})

test_that("plot.summary.matchit: runs for matched, unmatched, and subclass summaries", {
  local_null_device()

  m <- matchit(f_pl, data = lalonde)

  expect_silent(plot(summary(m)))
  expect_silent(plot(summary(m, un = FALSE)))
  expect_silent(plot(summary(m, interactions = TRUE)))
  expect_silent(plot(summary(matchit(f_pl, data = lalonde, method = NULL))))

  #Subclass summaries add one series of points per subclass
  ms <- matchit(f_pl, data = lalonde, method = "subclass", subclass = 3)
  expect_silent(plot(summary(ms)))
  expect_silent(plot(summary(ms, subclass = TRUE)))
})

test_that("plot.summary.matchit: needs standardized statistics", {
  local_null_device()

  s <- summary(matchit(f_pl, data = lalonde), standardize = FALSE)

  expect_err(plot(s), "not appropriate for unstandardized summary")
})

test_that("plot.summary.matchit: var.order, abs, threshold, and position", {
  local_null_device()

  m <- matchit(f_pl, data = lalonde)
  s <- summary(m)

  for (var.order in c("data", "matched", "unmatched", "alphabetical")) {
    expect_silent(plot(s, var.order = var.order))
  }

  expect_err(plot(s, var.order = "nonesuch"), "var.order")

  #Signed differences are plotted with thresholds on both sides
  expect_silent(plot(s, abs = FALSE))
  expect_silent(plot(s, threshold = NULL))
  expect_silent(plot(s, threshold = 0.1))
  expect_silent(plot(s, position = NULL))
  expect_silent(plot(s, position = "topright"))

  expect_err(plot(s, abs = NA), "abs")
})

test_that("plot.summary.matchit: var.order is rejected when its series is absent", {
  local_null_device()

  #`un = FALSE` drops the unmatched statistics, so they cannot be ordered by
  expect_err(plot(summary(matchit(f_pl, data = lalonde), un = FALSE),
                  var.order = "unmatched"),
             "`var.order` cannot be \"unmatched\" if `un = FALSE` in the call to `summary()`")

  #`method = NULL` produces no matched statistics
  expect_err(plot(summary(matchit(f_pl, data = lalonde, method = NULL)),
                  var.order = "matched"),
             "`var.order` cannot be \"matched\" if `method = NULL` in the original call")
})
