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
expect_matchit_snapshot <- function(m) {
  expect_snapshot_value(m$match.matrix, style = "json2")

  expect_snapshot_value(unname(round(m$weights, 8L)), style = "json2")

  subclass <- {
    if (is_null(m$subclass)) NULL
    else unname(as.integer(m$subclass))
  }

  expect_snapshot_value(subclass, style = "json2")

  invisible(m)
}

#Use regex to make strings invariant to white spaces
.w <- function(x) {
  gsub(" ", "(\\s+)", x, fixed = TRUE)
}

# ===== Condition helpers =====

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

#`matchit()` messages are built by cli, which hard-wraps them, so a regex spanning
#a wrap point silently fails to match. Use this on any expected message long enough
#to wrap. `expect_error(m, .w("..."))` is equivalent for short messages.
expect_matchit_condition <- function(expr, class = c("warning", "error"), pattern) {
  class <- match.arg(class)

  cnd <- tryCatch({
    withCallingHandlers(expr,
                        warning = function(w) {
                          if (class == "warning") {
                            stop(w)
                          }
                          invokeRestart("muffleWarning")
                        })
    NULL
  }, condition = function(c) c)

  expect(is_not_null(cnd),
         sprintf("No %s was signaled.", class),
         trace_env = rlang::caller_env())

  if (is_null(cnd)) {
    return(invisible(NULL))
  }

  #cli inserts newlines at its wrap width; collapse them before matching
  msg <- gsub("\\s+", " ", conditionMessage(cnd))

  expect_match(msg, pattern, fixed = TRUE)

  invisible(cnd)
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