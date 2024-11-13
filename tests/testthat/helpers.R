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
          #May not be duplicates incidentally; make sure examples induce duplicates
          expect_true(!isTRUE(all.equal(anyDuplicated(na.omit(m$match.matrix)), 0)))

          if (!is.null(attr(replace, "reuse.max"))) {
            expect_true(max(table(m$match.matrix)) <= attr(replace, "reuse.max"))
          }
        }
        else {
          expect_equal(anyDuplicated(na.omit(m$match.matrix)), 0)
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

  invisible(m)
}