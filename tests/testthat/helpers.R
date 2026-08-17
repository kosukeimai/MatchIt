#Helper functions for testing

expect_good_matchit <- function(m, expect_subclass = NULL, expect_distance = NULL,
                                expect_match.matrix = NULL, ratio = NULL, replace = NULL) {

  if (isTRUE(expect_subclass) && isTRUE(replace)) {
    stop("`expect_subclass` and `replace` cannot both be TRUE")
  }
  if (isTRUE(expect_match.matrix) && is.null(ratio)) {
    stop("`ratio` cannot be NULL when `expect_match.matrix` is TRUE")
  }

  expect_s3_class(m, "matchit")

  focal <- {
    if (is.null(m$focal)) switch(m$estimand, ATC = 0, 1)
    else m$focal
  }

  n <- length(m$treat)

  n1 <- sum(m$treat == focal)

  #Related to subclass
  if (!is.null(expect_subclass)) {
    if (expect_subclass) {
      expect_false(is.null(m$subclass))
      expect_length(m$subclass, n)
      expect_true(is.factor(m$subclass))
      expect_false(is.null(names(m$subclass)))
    }
    else {
      expect_null(m$subclass)
    }
  }

  #Related to match.matrix
  if (!is.null(expect_match.matrix)) {
    if (expect_match.matrix) {
      expect_false(is.null(m$match.matrix))
      expect_true(is.matrix(m$match.matrix))
      expect_true(is.character(m$match.matrix))
      expect_equal(nrow(m$match.matrix), n1)
      expect_equal(ncol(m$match.matrix), max(ratio[1], attr(ratio, "max.controls")))
      expect_false(is.null(rownames(m$match.matrix)))
      expect_false(any(rownames(m$match.matrix) %in% m$match.matrix))

      #Check no duplicates within each row
      expect_true(all(apply(m$match.matrix, 1, function(i) anyDuplicated(na.omit(i)) == 0)))

      if (!is.null(replace)) {
        if (replace) {
          #May not be duplicates incidentally; make sure examples induce duplicates;
          #if not expected to, set replace = FALSE
          expect_true(max(table(m$match.matrix)) > 0)

          if (!is.null(attr(replace, "reuse.max"))) {
            expect_true(max(table(m$match.matrix)) <= attr(replace, "reuse.max"))
          }
        }
        else {
          expect_true(max(table(m$match.matrix)) == 1)
        }
      }
    }
    else {
      expect_null(m$match.matrix)
    }
  }

  #Related to distance
  if (!is.null(expect_distance)) {
    if (expect_distance) {
      expect_false(is.null(m$distance))
      expect_length(m$distance, n)
      expect_true(is.numeric(m$distance))
      expect_false(anyNA(m$distance))
      expect_null(dim(m$distance))
      expect_false(is.null(names(m$distance)))
    }
    else {
      expect_null(m$distance)
    }
  }

  #Check that weights are equal when computed from subclasses or match.matrix
  if (!is.null(m$subclass) && !is.null(m$match.matrix)) {
    expect_equal(get_weights_from_mm(m$match.matrix, m$treat, focal),
                 get_weights_from_subclass(m$subclass, m$treat, m$estimand))
  }

  #Check that calipers work as expected
  if (!is.null(m$caliper)) {
    if (!is.null(m$match.matrix)) {
      expect_true(all(vapply(seq_along(m$caliper), function(i) {
        if (is.null(names(m$caliper))) {
          cc <- m$caliper[1L]
          xx <- m$distance
        }
        else if (names(m$caliper)[i] == "") {
          cc <- m$caliper[i]
          xx <- m$distance
        }
        else {
          cc <- m$caliper[i]
          xx <- setNames(m$X[[names(m$caliper)[i]]], names(m$treat))
        }

        all(vapply(rownames(m$match.matrix), function(t1) {
          all(vapply(na.omit(m$match.matrix[t1, ]), function(t0) {
            if (cc >= 0) {
              abs(xx[t1] - xx[t0]) <= cc
            }
            else {
              abs(xx[t1] - xx[t0]) > -cc
            }
          }, logical(1L)))
        }, logical(1L)))
      }, logical(1L))))
    }

    if (!is.null(m$subclass)) {
      expect_true(all(vapply(seq_along(m$caliper), function(i) {
        if (is.null(names(m$caliper))) {
          cc <- m$caliper[1L]
          xx <- m$distance
        }
        else if (names(m$caliper)[i] == "") {
          cc <- m$caliper[i]
          xx <- m$distance
        }
        else {
          cc <- m$caliper[i]
          xx <- setNames(m$X[[names(m$caliper)[i]]], names(m$treat))
        }

        all(vapply(which(m$treat == focal & !is.na(m$subclass)), function(t1) {
          all(vapply(which(m$treat != focal & !is.na(m$subclass) & m$subclass == m$subclass[t1]), function(t0) {
            if (cc >= 0) {
              abs(xx[t1] - xx[t0]) <= cc
            }
            else {
              abs(xx[t1] - xx[t0]) > -cc
            }
          }, logical(1L)))
        }, logical(1L)))
      }, logical(1L))))
    }
  }

  invisible(m)
}

#Pin the matched output of a `matchit()` call. `match.matrix` is recorded first so
#that snapshots taken before `weights` and `subclass` were added still match.
#`weights` needs its own pin because when `reuse.max > 1` there is no `subclass`
#and the weights are computed from `match.matrix` by a code path that no other
#expectation exercises; `subclass` needs one because several distinct subclass
#assignments can produce the same weights.
#
#These pins exist to hold results fixed across refactoring on one machine. They are
#not portable: the solver-backed methods record values that depend on the installed
#optmatch, quickmatch, or highs, and the propensity-score-based ones can shift on a
#different BLAS. Skipping on CRAN is deliberate. The skip lives here rather than at
#the top of each test so that the structural checks preceding it still run -- a
#failure recorded before a skip is still reported as a failure.
expect_matchit_snapshot <- function(m) {
  skip_on_cran()

  expect_snapshot_value(m$match.matrix, style = "json2")

  expect_snapshot_value(unname(round(m$weights, 8L)), style = "json2")

  subclass <- {
    if (is_null(m$subclass)) NULL
    else unname(as.integer(m$subclass))
  }

  expect_snapshot_value(subclass, style = "json2")

  invisible(m)
}

# ===== Condition helpers =====

#Messages from *arg* are built by cli, which hard-wraps them at the console width, so
#`conditionMessage()` contains newlines at positions that depend on that width and on
#the length of any interpolated values. A pattern spanning one of those breaks fails to
#match, silently and for a reason that has nothing to do with the behavior under test.
#This undoes the wrapping so an expected message can be written as one string and
#compared literally.
.collapse_cnd_message <- function(cnd) {
  gsub("\\s+", " ", trimws(conditionMessage(cnd)))
}

#Assert that `expr` signals a condition of `class`, and that each element of `pattern`
#appears in the message of one such condition. `pattern` may name several messages, so a
#call that warns more than once is asserted in a single expectation rather than by
#nesting `expect_warning()`s; omitting it asserts only that the condition fired.
#
#Matching is literal by default: *arg* messages are full of backticks, quotes, brackets,
#periods, and question marks, so treating them as regular expressions is almost always
#a mistake.
#
#Warnings and messages are caught with `withCallingHandlers()` so that `expr` runs to
#completion -- `tryCatch()` would unwind it, leaving anything assigned inside `expr`
#unset. Matching conditions are muffled; use `expect_no_unexpected_warning()` when the
#point is that nothing *else* warned.
expect_condition_message <- function(expr, class, pattern = NULL, fixed = TRUE) {
  matched <- list()
  err <- NULL

  handle <- function(cnd) {
    if (inherits(cnd, class)) {
      matched[[length(matched) + 1L]] <<- cnd
    }
  }

  tryCatch(
    withCallingHandlers(expr,
                        warning = function(w) {
                          handle(w)
                          invokeRestart("muffleWarning")
                        },
                        message = function(m) {
                          handle(m)
                          invokeRestart("muffleMessage")
                        }),
    error = function(e) {
      handle(e)

      if (!inherits(e, class)) {
        err <<- e
      }
    }
  )

  if (is_null(matched)) {
    #Report an unexpected error rather than only the absence of the expected condition;
    #otherwise a call that fails early looks identical to one that simply did not warn.
    expect(FALSE,
           sprintf("No condition of class %s was signaled.%s",
                   dQuote(class, FALSE),
                   if (is_null(err)) ""
                   else sprintf("\n  An error was signaled instead: %s",
                                encodeString(.collapse_cnd_message(err), quote = '"'))),
           trace_env = rlang::caller_env())

    return(invisible(NULL))
  }

  msgs <- vapply(matched, .collapse_cnd_message, character(1L))

  for (p in pattern) {
    expect(any(grepl(p, msgs, fixed = fixed)),
           sprintf("No %s message contained the expected text.\n  expected: %s\n  actual:   %s",
                   class, encodeString(p, quote = '"'),
                   toString(encodeString(msgs, quote = '"'))),
           trace_env = rlang::caller_env())
  }

  invisible(if (length(matched) == 1L) matched[[1L]] else matched)
}

#Counterparts to `arg::err()`, `arg::wrn()`, and `arg::msg()`, which signal rlang
#conditions. Matching on the rlang subclass rather than the base class also asserts
#that the condition came from *arg* (or another rlang-based signaller) rather than from
#base R, so a low-level failure that happens to mention the same words cannot pass for
#a proper input check.
expect_err <- function(expr, pattern = NULL, fixed = TRUE) {
  expect_condition_message(expr, "rlang_error", pattern, fixed)
}

expect_wrn <- function(expr, pattern = NULL, fixed = TRUE) {
  expect_condition_message(expr, "rlang_warning", pattern, fixed)
}

expect_msg <- function(expr, pattern = NULL, fixed = TRUE) {
  expect_condition_message(expr, "rlang_message", pattern, fixed)
}

#Assert that `expr` emits no warning other than ones matching `known`. Some
#specifications warn only for some samples or some solver versions, so they cannot
#simply be wrapped in `expect_warning()`; this accepts the known warnings and
#nothing else, so an unrelated condition still fails the test, and it keeps working
#unchanged if the known warning stops firing.
expect_no_unexpected_warning <- function(expr, known = character()) {
  ws <- character()

  val <- withCallingHandlers(expr,
                             warning = function(w) {
                               ws <<- c(ws, conditionMessage(w))
                               invokeRestart("muffleWarning")
                             })

  unexpected <- {
    if (is_null(known)) ws
    else ws[!vapply(ws, function(w) any(vapply(known, grepl, logical(1L),
                                               x = w, fixed = TRUE)),
                    logical(1L))]
  }

  expect(is_null(unexpected),
         sprintf("Unexpected %s signaled: %s%s",
                 ngettext(length(unexpected), "warning was", "warnings were"),
                 toString(dQuote(gsub("\\s+", " ", unexpected), FALSE)),
                 if (is_null(known)) ""
                 else sprintf(" (allowed: %s)", toString(dQuote(known, FALSE)))),
         trace_env = rlang::caller_env())

  invisible(val)
}

# ===== Comparison helpers =====

#`expect_equal()`'s complement, for asserting that something actually changed.
expect_not_equal <- function(object, expected, ...,
                             tolerance = if (edition_get() >= 3) testthat_tolerance(),
                             info = NULL, label = NULL, expected.label = NULL) {

  act <- quasi_label(rlang::enquo(object), label, arg = "object")
  exp <- quasi_label(rlang::enquo(expected), expected.label, arg = "expected")

  comp <- {
    if (is_null(tolerance)) waldo::compare(act$val, exp$val, ..., x_arg = "actual",
                                           y_arg = "expected")
    else waldo::compare(act$val, exp$val, ..., tolerance = tolerance,
                        x_arg = "actual", y_arg = "expected")
  }

  expect(length(comp) > 0L,
         sprintf("%s (`actual`) is equal to %s (`expected`).", act$lab, exp$lab),
         info = info, trace_env = rlang::caller_env())

  invisible(act$val)
}

#Matching is supposed to improve balance; this is the weakest check that it did
#anything useful, and it catches sign errors and mixed-up weights that structural
#checks and snapshots both pass over.
expect_balance_improved <- function(m, ...) {
  skip_if_not_installed("cobalt")

  s.weights <- m$s.weights %or% rep.int(1, length(m$treat))

  matched <- abs(cobalt::col_w_smd(m$X, m$treat, m$weights, s.weights = s.weights))
  unmatched <- abs(cobalt::col_w_smd(m$X, m$treat, s.weights = s.weights))

  expect_lt(mean(matched, na.rm = TRUE), mean(unmatched, na.rm = TRUE), ...)

  invisible(m)
}

# ===== Fixtures =====

#Set a fixed proportion of each named column to NA without disturbing the RNG
#stream of the calling test.
inject_missingness <- function(data, cols, prop = 0.1, seed = 4321) {
  old_seed <- {
    if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv())
    else NULL
  }

  set.seed(seed)

  for (col in cols) {
    is.na(data[[col]]) <- sample(nrow(data), round(prop * nrow(data)))
  }

  if (is_null(old_seed)) {
    rm(".Random.seed", envir = globalenv())
  }
  else {
    assign(".Random.seed", old_seed, envir = globalenv())
  }

  data
}
