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
#
#Taken from cobalt's test helpers, which solve this problem already.
#
#Messages raised by `arg::err()`/`arg::wrn()`/`arg::msg()` are formatted by cli,
#which capitalizes the first letter, appends a period, converts inline markup
#(e.g., `{.arg x}` to `` `x` ``), and hard-wraps the result to the console width.
#The wrapping means a literal `fixed = TRUE` match against a long message fails,
#and building a regex from the message requires escaping every metacharacter that
#cli may have introduced (`[`, `]`, `{`, `}`, `+`, `?`, `*`, `|`).
#
#`expect_err()`, `expect_wrn()`, and `expect_msg()` sidestep both problems: they
#collapse all whitespace in the *observed* message and then do a literal substring
#match. Never copy a source string from `R/` into these -- use the rendered text.

#Collapse runs of whitespace so matching is invariant to how cli wrapped the message.
squish <- function(x) {
  gsub("\\s+", " ", trimws(paste(x, collapse = " ")))
}

#Shared back end for the three expectations below.
.expect_cnd_text <- function(cnd, text, what) {
  if (is.null(cnd)) {
    testthat::fail(sprintf("Expected %s, but none was signaled.", what))
    return(invisible(NULL))
  }

  if (is.null(text)) {
    testthat::succeed()
    return(invisible(cnd))
  }

  msg <- squish(conditionMessage(cnd))

  #Matching is case-insensitive because cli capitalizes the first letter of every
  #message, which would otherwise make any substring starting at the beginning of
  #the message fail. The point is to identify which message was raised, not to
  #pin its capitalization.
  testthat::expect_true(
    grepl(tolower(text), tolower(msg), fixed = TRUE),
    info = sprintf("%s message did not contain the expected text.\n  expected: %s\n  actual:   %s",
                   what, encodeString(text, quote = "\""), encodeString(msg, quote = "\""))
  )

  invisible(cnd)
}

#Each of the three expectations below asserts that one condition was signaled.
#Anything else the expression prints or signals is incidental to that assertion, and
#letting it through only buries the reporter's own output, so the other two condition
#classes are muffled and printed output is discarded.
#
#`capture.output()` evaluates its argument in the calling frame, so the promise is
#still forced there and assignments inside the expression still take effect there.

#Expect an error whose message contains `text` (literal, whitespace-insensitive).
expect_err <- function(object, text = NULL) {
  cnd <- NULL

  utils::capture.output(
    withCallingHandlers(cnd <- tryCatch({
      force(object)
      NULL
    }, error = function(e) e),
    message = function(m) invokeRestart("muffleMessage"),
    warning = function(w) invokeRestart("muffleWarning")))

  .expect_cnd_text(cnd, text, "error")
}

#Expect a warning whose message contains `text`. The expression still runs to
#completion, so assignments inside it take effect in the calling environment.
expect_wrn <- function(object, text = NULL) {
  cnd <- NULL

  utils::capture.output(
    withCallingHandlers(force(object),
                        message = function(m) invokeRestart("muffleMessage"),
                        warning = function(w) {
                          if (is.null(cnd)) cnd <<- w
                          invokeRestart("muffleWarning")
                        }))

  .expect_cnd_text(cnd, text, "warning")
}

#Expect a message whose message contains `text`. As with `expect_wrn()`, the
#expression runs to completion.
expect_msg <- function(object, text = NULL) {
  cnd <- NULL

  utils::capture.output(
    withCallingHandlers(force(object),
                        warning = function(w) invokeRestart("muffleWarning"),
                        message = function(m) {
                          if (is.null(cnd)) cnd <<- m
                          invokeRestart("muffleMessage")
                        }))

  .expect_cnd_text(cnd, text, "message")
}

#Assert that `expr` emits no warning other than ones matching `known`. Some
#specifications warn only for some samples or some solver versions, so they cannot
#simply be wrapped in `expect_wrn()`; this accepts the known warnings and nothing
#else, so an unrelated condition still fails the test, and it keeps working
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
                 toString(dQuote(squish_each(unexpected), FALSE)),
                 if (is_null(known)) ""
                 else sprintf(" (allowed: %s)", toString(dQuote(known, FALSE)))),
         trace_env = rlang::caller_env())

  invisible(val)
}

#`squish()` collapses a vector into one string; this collapses each element.
squish_each <- function(x) {
  vapply(x, squish, character(1L), USE.NAMES = FALSE)
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

#Send plots to a null device for the duration of the calling test, so the plotting
#functions can be exercised without opening a window or writing a file. Taken from
#cobalt's test helpers.
local_null_device <- function(.env = parent.frame()) {
  grDevices::pdf(NULL)

  #`rlang::defer()` is not exported, so register the cleanup as an `on.exit()`
  #expression in the calling frame instead. This needs no extra dependency.
  do.call(base::on.exit,
          list(quote(grDevices::dev.off()), add = TRUE, after = FALSE),
          envir = .env)

  invisible(NULL)
}

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
