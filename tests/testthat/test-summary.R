# Tests for summary.matchit() and summary.matchit.subclass(), plus their print methods.
#
# The balance matrices are pinned rather than the printed output: printed tables depend
# on `getOption("digits")` and `OutDec` and on the column-formatting helper, so a text
# snapshot would break for reasons unrelated to the statistics. The printing itself is
# checked structurally, by asserting which blocks appear.

data("lalonde", package = "MatchIt")

f_sum <- treat ~ age + educ + race + married + nodegree + re74 + re75

lalonde_sw <- seq(0.5, 2, length.out = nrow(lalonde))

#Pin the numbers a `summary()` reports. Rounded for the same reason the matching
#snapshots are: the last bits are platform noise, not behavior.
expect_summary_snapshot <- function(s) {
  skip_on_cran()

  for (component in c("sum.all", "sum.matched", "sum.across", "reduction")) {
    if (is_not_null(s[[component]])) {
      expect_snapshot_value(round(unclass(s[[component]]), 6L), style = "json2")
    }
  }

  expect_snapshot_value(round(unclass(s$nn), 6L), style = "json2")

  invisible(s)
}

# ===== structure =====

test_that("summary: components and class", {
  m <- matchit(f_sum, data = lalonde, method = "nearest")
  s <- summary(m)

  expect_s3_class(s, "summary.matchit")
  expect_named(s, c("call", "nn", "sum.all", "sum.matched", "reduction"))

  expect_true(is.matrix(s$sum.all))
  expect_true(is.matrix(s$sum.matched))
  expect_identical(rownames(s$sum.all), rownames(s$sum.matched))

  expect_identical(colnames(s$sum.matched),
                   c("Means Treated", "Means Control", "Std. Mean Diff.",
                     "Var. Ratio", "eCDF Mean", "eCDF Max", "Std. Pair Dist."))

  expect_summary_snapshot(s)
})

test_that("summary: the nn table counts what it says", {
  m <- matchit(f_sum, data = lalonde, method = "nearest")
  nn <- summary(m)$nn

  expect_identical(rownames(nn),
                   c("All (ESS)", "All", "Matched (ESS)", "Matched",
                     "Unmatched", "Discarded"))
  expect_identical(colnames(nn), c("Control", "Treated"))

  for (g in c(Control = 0L, Treated = 1L)) {
    col <- names(which(c(Control = 0L, Treated = 1L) == g))

    expect_equal(nn["All", col], sum(m$treat == g), ignore_attr = TRUE)
    expect_equal(nn["Matched", col], sum(m$weights > 0 & m$treat == g),
                 ignore_attr = TRUE)
    expect_equal(nn["Unmatched", col],
                 sum(m$weights == 0 & !m$discarded & m$treat == g),
                 ignore_attr = TRUE)
    expect_equal(nn["Discarded", col], sum(m$discarded & m$treat == g),
                 ignore_attr = TRUE)
  }

  #With no sampling weights and 1:1 matching the ESS rows equal the raw counts
  expect_equal(nn["All (ESS)", ], nn["All", ])
  expect_equal(nn["Matched (ESS)", ], nn["Matched", ])
})

test_that("summary: discarded units are counted separately from unmatched", {
  m <- matchit(f_sum, data = lalonde, method = "nearest", discard = "both")
  nn <- summary(m)$nn

  expect_gt(sum(nn["Discarded", ]), 0L)
  expect_equal(sum(nn["Discarded", ]), sum(m$discarded), ignore_attr = TRUE)

  #Every unit is accounted for exactly once
  expect_equal(sum(nn[c("Matched", "Unmatched", "Discarded"), ]),
               nrow(lalonde), ignore_attr = TRUE)
})

test_that("summary: s.weights make the ESS rows differ from the counts", {
  m <- matchit(f_sum, data = lalonde, method = "nearest", s.weights = lalonde_sw)
  nn <- summary(m)$nn

  expect_not_equal(nn["All (ESS)", ], nn["All", ])
  expect_lt(nn["All (ESS)", "Control"], nn["All", "Control"])

  #The raw counts are unaffected by weighting
  m0 <- matchit(f_sum, data = lalonde, method = "nearest")
  expect_equal(nn["All", ], summary(m0)$nn["All", ])

  expect_summary_snapshot(summary(m))
})

# ===== options =====

test_that("summary: un = FALSE drops the unmatched-sample block", {
  m <- matchit(f_sum, data = lalonde, method = "nearest")

  expect_null(summary(m, un = FALSE)$sum.all)
  expect_false(is_null(summary(m, un = TRUE)$sum.all))

  #...and the corresponding block disappears from the printout
  expect_true(any(grepl("All Data", capture.output(print(summary(m))), fixed = TRUE)))
  expect_false(any(grepl("All Data", capture.output(print(summary(m, un = FALSE))),
                         fixed = TRUE)))
})

test_that("summary: un = FALSE is overridden when there was no matching", {
  #With `method = NULL` there is nothing to compare against, so the all-data block is
  #reported whatever `un` says.
  m <- matchit(f_sum, data = lalonde, method = NULL)

  expect_false(is_null(summary(m, un = FALSE)$sum.all))
  expect_null(summary(m)$sum.matched)
})

test_that("summary: improvement controls whether reduction is computed", {
  m <- matchit(f_sum, data = lalonde, method = "nearest")

  expect_null(summary(m)$reduction)

  red <- summary(m, improvement = TRUE)$reduction
  expect_true(is.matrix(red))
  expect_identical(nrow(red), nrow(summary(m)$sum.all))
  expect_identical(colnames(red),
                   c("Std. Mean Diff.", "Var. Ratio", "eCDF Mean", "eCDF Max",
                     "Std. Pair Dist."))

  expect_true(any(grepl("Improvement",
                        capture.output(print(summary(m, improvement = TRUE))))))
  expect_false(any(grepl("Improvement", capture.output(print(summary(m))))))
})

test_that("summary: standardize changes both the column names and the values", {
  m <- matchit(f_sum, data = lalonde, method = "nearest")

  s_std <- summary(m, standardize = TRUE)$sum.matched
  s_raw <- summary(m, standardize = FALSE)$sum.matched

  expect_identical(colnames(s_raw),
                   c("Means Treated", "Means Control", "Mean Diff.", "Var. Ratio",
                     "eQQ Mean", "eQQ Max", "Pair Dist."))

  #The group means are the same; the difference columns are not
  expect_equal(s_std[, 1:2], s_raw[, 1:2], ignore_attr = TRUE)
  expect_not_equal(unname(s_std[, 3L]), unname(s_raw[, 3L]))

  expect_summary_snapshot(summary(m, standardize = FALSE))
})

test_that("summary: pair.dist = FALSE blanks the pair-distance column only", {
  m <- matchit(f_sum, data = lalonde, method = "nearest")

  with_pd <- summary(m, pair.dist = TRUE)$sum.matched
  without <- summary(m, pair.dist = FALSE)$sum.matched

  expect_false(all(is.na(with_pd[, "Std. Pair Dist."])))
  expect_true(all(is.na(without[, "Std. Pair Dist."])))
  expect_equal(with_pd[, 1:6], without[, 1:6])
})

test_that("summary: interactions adds squares and interactions", {
  m <- matchit(f_sum, data = lalonde, method = "nearest")

  s0 <- summary(m)
  s1 <- summary(m, interactions = TRUE)

  expect_gt(nrow(s1$sum.all), nrow(s0$sum.all))

  #The original rows are still present and unchanged
  expect_equal(s1$sum.all[rownames(s0$sum.all), ], s0$sum.all)

  expect_summary_snapshot(s1)
})

test_that("summary: addlvariables accepts a formula, a character vector, and a data frame", {
  m <- matchit(f_sum, data = lalonde, method = "nearest")
  n0 <- nrow(summary(m)$sum.all)

  #A formula, evaluated in the data or the environment
  s_f <- summary(m, addlvariables = ~ I(age^2))
  expect_identical(nrow(s_f$sum.all), n0 + 1L)

  #A data frame with one row per original unit
  s_d <- summary(m, addlvariables = data.frame(twice_age = lalonde$age * 2))
  expect_identical(nrow(s_d$sum.all), n0 + 1L)
  expect_true("twice_age" %in% rownames(s_d$sum.all))

  #Names of variables already in the model add nothing new
  expect_identical(nrow(summary(m, addlvariables = c("age", "educ"))$sum.all), n0)

  expect_summary_snapshot(s_f)
})

test_that("summary: addlvariables can name variables found only in `data`", {
  m <- matchit(f_sum, data = lalonde, method = "nearest")

  lalonde_extra <- lalonde
  lalonde_extra$extra <- lalonde$age * 3

  s <- summary(m, addlvariables = ~ extra, data = lalonde_extra)
  expect_true("extra" %in% rownames(s$sum.all))

  #Without `data` the variable cannot be found
  expect_error(summary(m, addlvariables = ~ extra))
})

test_that("summary: addlvariables accepts a matrix", {
  m <- matchit(f_sum, data = lalonde, method = "nearest")
  n0 <- nrow(summary(m)$sum.all)

  #A numeric matrix is coerced to a data frame; columns already in the model add
  #nothing new, and unnamed columns are named as `as.data.frame()` names them
  s <- summary(m, addlvariables = cbind(age = lalonde$age,
                                        twice_age = lalonde$age * 2))
  expect_identical(nrow(s$sum.all), n0 + 1L)
  expect_true("twice_age" %in% rownames(s$sum.all))

  s_unnamed <- summary(m, addlvariables = matrix(lalonde$age * 3, ncol = 1L))
  expect_true("V1" %in% rownames(s_unnamed$sum.all))

  #A wrongly sized matrix is an error naming `addlvariables`, not `data`
  expect_err(summary(m, addlvariables = as.matrix(lalonde[1:5, "age", drop = FALSE])),
             "variables specified in `addlvariables` must have the same number")
})

test_that("summary: addlvariables in an unaccepted form is an error naming addlvariables", {
  m <- matchit(f_sum, data = lalonde, method = "nearest")

  expect_err(summary(m, addlvariables = as.list(lalonde)),
             "the argument to `addlvariables` must be in one of the accepted forms")
  expect_err(summary(m, addlvariables = 1:614),
             "the argument to `addlvariables` must be in one of the accepted forms")
})

# ===== subclassification =====

test_that("summary: subclassification has its own class and components", {
  m <- matchit(f_sum, data = lalonde, method = "subclass", subclass = 4)
  s <- summary(m)

  expect_s3_class(s, "summary.matchit.subclass")
  expect_s3_class(s, "summary.matchit")
  expect_named(s, c("call", "sum.all", "sum.across", "sum.subclass", "reduction",
                    "qn", "nn"))

  #`qn` counts units per subclass per treatment group, plus a total column
  expect_identical(nrow(s$qn), 3L)
  expect_identical(ncol(s$qn), nlevels(m$subclass) + 1L)

  expect_summary_snapshot(s)
})

test_that("summary: the subclass argument selects which subclasses to report", {
  m <- matchit(f_sum, data = lalonde, method = "subclass", subclass = 4)

  expect_length(summary(m, subclass = FALSE)$sum.subclass, 0L)
  expect_length(summary(m, subclass = TRUE)$sum.subclass, 4L)

  expect_named(summary(m, subclass = TRUE)$sum.subclass,
               paste("Subclass", 1:4))

  #A subset of indices, in the order given
  s2 <- summary(m, subclass = 2)
  expect_length(s2$sum.subclass, 1L)
  expect_named(s2$sum.subclass, "Subclass 2")
  expect_equal(s2$sum.subclass[[1L]],
               summary(m, subclass = TRUE)$sum.subclass[["Subclass 2"]])

  expect_length(summary(m, subclass = c(1, 3))$sum.subclass, 2L)
})

test_that("summary: an out-of-range subclass index is an error", {
  m <- matchit(f_sum, data = lalonde, method = "subclass", subclass = 4)

  expect_err(summary(m, subclass = 99),
             "`subclass` should be TRUE, FALSE, or a vector of subclass indices")
})

# ===== across methods =====

test_that("summary: works for every matching method", {
  fits <- list(
    nearest = matchit(f_sum, data = lalonde, method = "nearest"),
    `nearest replace` = matchit(f_sum, data = lalonde, method = "nearest",
                                replace = TRUE),
    cem = matchit(f_sum, data = lalonde, method = "cem"),
    exact = matchit(treat ~ age + educ + race, data = lalonde, method = "exact"),
    cardinality = matchit(f_sum, data = lalonde, method = "cardinality"),
    `no matching` = matchit(f_sum, data = lalonde, method = NULL)
  )

  if (rlang::is_installed("optmatch")) {
    fits$full <- matchit(f_sum, data = lalonde, method = "full")
    fits$optimal <- matchit(f_sum, data = lalonde, method = "optimal")
  }

  for (nm in names(fits)) {
    s <- summary(fits[[nm]])
    expect_s3_class(s, "summary.matchit")
    expect_true(is.matrix(s$nn))
    expect_no_error(invisible(capture.output(print(s))))
  }

  expect_summary_snapshot(summary(fits$cem))
})

test_that("summary: pair distances come from the strata, not from match.matrix", {
  #`Std. Pair Dist.` is computed within subclasses, so it is reported for methods that
  #produce a `subclass` even when they produce no `match.matrix`, and is absent only
  #when there is neither.
  m_cem <- matchit(f_sum, data = lalonde, method = "cem")
  expect_null(m_cem$match.matrix)
  expect_false(all(is.na(summary(m_cem)$sum.matched[, "Std. Pair Dist."])))

  #Cardinality matching selects a subset without forming strata or pairs
  m_card <- matchit(f_sum, data = lalonde, method = "cardinality")
  expect_null(m_card$subclass)
  expect_null(m_card$match.matrix)
  expect_true(all(is.na(summary(m_card)$sum.matched[, "Std. Pair Dist."])))

  #Under exact matching every unit in a stratum has identical covariate values, so
  #every within-stratum distance is zero
  m_exact <- matchit(treat ~ age + educ + race, data = lalonde, method = "exact")
  pd <- summary(m_exact)$sum.matched[, "Std. Pair Dist."]
  expect_equal(unname(pd[!is.na(pd)]), rep(0, sum(!is.na(pd))))
})

# ===== printing =====

test_that("summary: print returns its input invisibly and respects digits", {
  m <- matchit(f_sum, data = lalonde, method = "nearest")
  s <- summary(m)

  expect_no_error(invisible(capture.output(print(s, digits = 2))))
  expect_not_equal(capture.output(print(s, digits = 2)),
                   capture.output(print(s, digits = 5)))

  expect_err(print(s, digits = "two"))
})

test_that("summary: print shows the call and the sample sizes", {
  m <- matchit(f_sum, data = lalonde, method = "nearest")
  out <- capture.output(print(summary(m)))

  expect_true(any(grepl("Call:", out, fixed = TRUE)))
  expect_true(any(grepl("Sample Sizes", out, fixed = TRUE)))
  expect_true(any(grepl("Matched Data", out, fixed = TRUE)))
})

test_that("summary: subclass print shows the per-subclass blocks only when asked", {
  m <- matchit(f_sum, data = lalonde, method = "subclass", subclass = 4)

  out_no <- capture.output(print(summary(m, subclass = FALSE)))
  out_yes <- capture.output(print(summary(m, subclass = TRUE)))

  expect_false(any(grepl("Subclass 1", out_no, fixed = TRUE)))
  expect_true(any(grepl("Subclass 1", out_yes, fixed = TRUE)))
})

test_that("summary: no unexpected conditions", {
  m <- matchit(f_sum, data = lalonde, method = "nearest")

  expect_no_unexpected_warning(summary(m))
  expect_no_unexpected_warning(summary(m, interactions = TRUE, improvement = TRUE))
})
