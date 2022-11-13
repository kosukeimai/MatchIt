#' Add sampling weights to a `matchit` object
#'
#' @description
#' Adds sampling weights to a `matchit` object so that they are
#' incorporated into balance assessment and creation of the weights. This would
#' typically only be used when an argument to `s.weights` was not supplied
#' to [matchit()] (i.e., because they were not to be included in the estimation
#' of the propensity score) but sampling weights are required for generalizing
#' an effect to the correct population. Without adding sampling weights to the
#' `matchit` object, balance assessment tools (i.e., [summary.matchit()]
#' and [plot.matchit()]) will not calculate balance statistics correctly, and
#' the weights produced by [match.data()] and [get_matches()] will not
#' incorporate the sampling weights.
#'
#' @param m a `matchit` object; the output of a call to [matchit()],
#' typically with the `s.weights` argument unspecified.
#' @param s.weights an numeric vector of sampling weights to be added to the
#' `matchit` object. Can also be specified as a string containing the name
#' of variable in `data` to be used or a one-sided formula with the
#' variable on the right-hand side (e.g., `~ SW`).
#' @param data a data frame containing the sampling weights if given as a
#' string or formula. If unspecified, `add_s.weights()` will attempt to find
#' the dataset using the environment of the `matchit` object.
#'
#' @return a `matchit` object with an `s.weights` component
#' containing the supplied sampling weights. If `s.weights = NULL`, the original
#' `matchit` object is returned.
#'
#' @author Noah Greifer
#'
#' @seealso [matchit()]; [match.data()]
#'
#' @examples
#'
#' data("lalonde")
#'
#' # Generate random sampling weights, just
#' # for this example
#' sw <- rchisq(nrow(lalonde), 2)
#'
#' # NN PS match using logistic regression PS that doesn't
#' # include sampling weights
#' m.out <- matchit(treat ~ age + educ + race + nodegree +
#'                    married  + re74 + re75, data = lalonde)
#'
#' m.out
#'
#' # Add s.weights to the matchit object
#' m.out <- add_s.weights(m.out, sw)
#'
#' m.out #note additional output
#'
#' # Check balance; note that sample sizes incorporate
#' # s.weights
#' summary(m.out, improvement = FALSE)
#'
#' @export
add_s.weights <- function(m,
                          s.weights = NULL,
                          data = NULL) {

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
