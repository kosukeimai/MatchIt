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

  n <- length(m$treat)
  n1 <- sum(m$treat == 1)

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
    expect_equal(get_weights_from_mm(m$match.matrix, m$treat, switch(m$estimand, ATC = 0, 1)),
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

        all(vapply(which(m$treat == 1 & !is.na(m$subclass)), function(t1) {
          all(vapply(which(m$treat == 0 & !is.na(m$subclass) & m$subclass == m$subclass[t1]), function(t0) {
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