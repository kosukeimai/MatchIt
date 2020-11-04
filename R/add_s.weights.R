add_s.weights <- function(m, s.weights = NULL, data = NULL) {

  if (!inherits(m, "matchit")) stop("'m' must be a matchit object, the output of a call to matchit().", call. = FALSE)

  if (!is.null(s.weights)) {
    if (!is.numeric(s.weights)) {
      if (is.null(data)) {
        if (!is.null(m$model)) {
          env <- attributes(terms(m$model))$.Environment
        } else {
          env <- parent.frame()
        }
        data <- eval(m$call$data, envir = env)
        if (length(data) == 0) stop("A dataset could not be found. Please supply an argument to 'data' containing the original dataset used in the matching.", call. = FALSE)
      }
      else {
        if (!is.data.frame(data)) {
          if (is.matrix(data)) data <- as.data.frame.matrix(data)
          else stop("'data' must be a data frame.", call. = FALSE)
        }
        if (nrow(data) != length(m$treat)) {
          stop("'data' must have as many rows as there were units in the original call to matchit().", call. = FALSE)
        }
      }

      if (is.character(s.weights)) {
        if (is.null(data) || !is.data.frame(data)) {
          stop("If 's.weights' is specified a string, a data frame containing the named variable must be supplied to 'data'.", call. = FALSE)
        }
        if (!all(s.weights %in% names(data))) {
          stop("The name supplied to 's.weights' must be a variable in 'data'.", call. = FALSE)
        }
        s.weights.form <- reformulate(s.weights)
        s.weights <- model.frame(s.weights.form, data, na.action = "na.pass")
        if (ncol(s.weights) != 1) stop("'s.weights' can only contain one named variable.", call. = FALSE)
        s.weights <- s.weights[[1]]
      }
      else if (inherits(s.weights, "formula")) {
        s.weights.form <- update(s.weights, NULL ~ .)
        s.weights <- model.frame(s.weights.form, data, na.action = "na.pass")
        if (ncol(s.weights) != 1) stop("'s.weights' can only contain one named variable.", call. = FALSE)
        s.weights <- s.weights[[1]]
      }
      else {
        stop("'s.weights' must be supplied as a numeric vector, string, or one-sided formula.", call. = FALSE)
      }
    }

    if (anyNA(s.weights)) stop("Missing values are not allowed in 's.weights'.", call. = FALSE)
    if (length(s.weights) != length(m$treat)) stop("'s.weights' must be the same length as the treatment vector.", call. = FALSE)

    names(s.weights) <- names(m$treat)

    attr(s.weights, "in_ps") <- isTRUE(all.equal(s.weights, m$s.weights))

    m$s.weights <- s.weights

    m$nn <- nn(m$treat, m$weights, m$discarded, s.weights)
  }

  return(m)
}