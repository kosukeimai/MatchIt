#' Propensity scores and other distance measures
#' @name distance
#' @aliases distance
#' @usage NULL
#'
#' @description
#' Several matching methods require or can involve the distance between treated
#' and control units. Options include the Mahalanobis distance, propensity
#' score distance, or distance between user-supplied values. Propensity scores
#' are also used for common support via the `discard` options and for
#' defining calipers. This page documents the options that can be supplied to
#' the `distance` argument to [matchit()].
#'
#' @section Allowable options:
#'
#' There are four ways to specify the `distance` argument: 1) as a string containing the name of a method for
#' estimating propensity scores, 2) as a string containing the name of a method
#' for computing pairwise distances from the covariates, 3) as a vector of
#' values whose pairwise differences define the distance between units, or 4)
#' as a distance matrix containing all pairwise distances. The options are
#' detailed below.
#'
#' ## Propensity score estimation methods
#'
#' When `distance` is specified as the name of a method for estimating propensity scores
#' (described below), a propensity score is estimated using the variables in
#' `formula` and the method corresponding to the given argument. This
#' propensity score can be used to compute the distance between units as the
#' absolute difference between the propensity scores of pairs of units.
#' Propensity scores can also be used to create calipers and common support
#' restrictions, whether or not they are used in the actual distance measure
#' used in the matching, if any.
#'
#' In addition to the `distance` argument, two other arguments can be
#' specified that relate to the estimation and manipulation of the propensity
#' scores. The `link` argument allows for different links to be used in
#' models that require them such as generalized linear models, for which the
#' logit and probit links are allowed, among others. In addition to specifying
#' the link, the `link` argument can be used to specify whether the
#' propensity score or the linearized version of the propensity score should be
#' used; by specifying `link = "linear.{link}"`, the linearized version
#' will be used.
#'
#' The `distance.options` argument can also be specified, which should be
#' a list of values passed to the propensity score-estimating function, for
#' example, to choose specific options or tuning parameters for the estimation
#' method. If `formula`, `data`, or `verbose` are not supplied
#' to `distance.options`, the corresponding arguments from
#' `matchit()` will be automatically supplied. See the Examples for
#' demonstrations of the uses of `link` and `distance.options`. When
#' `s.weights` is supplied in the call to `matchit()`, it will
#' automatically be passed to the propensity score-estimating function as the
#' `weights` argument unless otherwise described below.
#'
#' The following methods for estimating propensity scores are allowed:
#'
#' \describe{
#' \item{`"glm"`}{ The propensity scores are estimated using
#' a generalized linear model (e.g., logistic regression). The `formula`
#' supplied to `matchit()` is passed directly to [glm()], and
#' [predict.glm()] is used to compute the propensity scores. The `link`
#' argument can be specified as a link function supplied to [binomial()], e.g.,
#' `"logit"`, which is the default. When `link` is prepended by
#' `"linear."`, the linear predictor is used instead of the predicted
#' probabilities. `distance = "glm"` with `link = "logit"` (logistic
#' regression) is the default in `matchit()`. (This used to be able to be requested as `distance = "ps"`, which still works.)}
#' \item{`"gam"`}{
#' The propensity scores are estimated using a generalized additive model. The
#' `formula` supplied to `matchit()` is passed directly to
#' \pkgfun{mgcv}{gam}, and \pkgfun{mgcv}{predict.gam} is used to compute the propensity
#' scores. The `link` argument can be specified as a link function
#' supplied to [binomial()], e.g., `"logit"`, which is the default. When
#' `link` is prepended by `"linear."`, the linear predictor is used
#' instead of the predicted probabilities. Note that unless the smoothing
#' functions \pkgfun{mgcv}{s}, \pkgfun{mgcv}{te}, \pkgfun{mgcv}{ti}, or \pkgfun{mgcv}{t2} are
#' used in `formula`, a generalized additive model is identical to a
#' generalized linear model and will estimate the same propensity scores as
#' `glm()`. See the documentation for \pkgfun{mgcv}{gam},
#' \pkgfun{mgcv}{formula.gam}, and \pkgfun{mgcv}{gam.models} for more information on
#' how to specify these models. Also note that the formula returned in the
#' `matchit()` output object will be a simplified version of the supplied
#' formula with smoothing terms removed (but all named variables present). }
#' \item{`"gbm"`}{ The propensity scores are estimated using a
#' generalized boosted model. The `formula` supplied to `matchit()`
#' is passed directly to \pkgfun{gbm}{gbm}, and \pkgfun{gbm}{predict.gbm} is used to
#' compute the propensity scores. The optimal tree is chosen using 5-fold
#' cross-validation by default, and this can be changed by supplying an
#' argument to `method` to `distance.options`; see \pkgfun{gbm}{gbm.perf}
#' for details. The `link` argument can be specified as `"linear"` to
#' use the linear predictor instead of the predicted probabilities. No other
#' links are allowed. The tuning parameter defaults differ from
#' `gbm::gbm()`; they are as follows: `n.trees = 1e4`,
#' `interaction.depth = 3`, `shrinkage = .01`, `bag.fraction = 1`, `cv.folds = 5`, `keep.data = FALSE`. These are the same
#' defaults as used in *WeightIt* and *twang*, except for
#' `cv.folds` and `keep.data`. Note this is not the same use of
#' generalized boosted modeling as in *twang*; here, the number of trees is
#' chosen based on cross-validation or out-of-bag error, rather than based on
#' optimizing balance. \pkg{twang} should not be cited when using this method
#' to estimate propensity scores. Note that because there is a random component to choosing the tuning
#' parameter, results will vary across runs unless a [seed][set.seed] is
#' set.}
#' \item{`"lasso"`, `"ridge"`, `"elasticnet"`}{
#'   The propensity
#' scores are estimated using a lasso, ridge, or elastic net model,
#' respectively. The `formula` supplied to `matchit()` is processed
#' with [model.matrix()] and passed to \pkgfun{glmnet}{cv.glmnet}, and
#' \pkgfun{glmnet}{predict.cv.glmnet} is used to compute the propensity scores. The
#' `link` argument can be specified as a link function supplied to
#' [binomial()], e.g., `"logit"`, which is the default. When `link`
#' is prepended by `"linear."`, the linear predictor is used instead of
#' the predicted probabilities. When `link = "log"`, a Poisson model is
#' used. For `distance = "elasticnet"`, the `alpha` argument, which
#' controls how to prioritize the lasso and ridge penalties in the elastic net,
#' is set to .5 by default and can be changed by supplying an argument to
#' `alpha` in `distance.options`. For `"lasso"` and
#' `"ridge"`, `alpha` is set to 1 and 0, respectively, and cannot be
#' changed. The `cv.glmnet()` defaults are used to select the tuning
#' parameters and generate predictions and can be modified using
#' `distance.options`. If the `s` argument is passed to
#' `distance.options`, it will be passed to `predict.cv.glmnet()`.
#' Note that because there is a random component to choosing the tuning
#' parameter, results will vary across runs unless a [seed][set.seed] is
#' set. }
#' \item{`"rpart"`}{ The propensity scores are estimated using a
#' classification tree. The `formula` supplied to `matchit()` is
#' passed directly to \pkgfun{rpart}{rpart}, and \pkgfun{rpart}{predict.rpart} is used
#' to compute the propensity scores. The `link` argument is ignored, and
#' predicted probabilities are always returned as the distance measure. }
#' \item{`"randomforest"`}{ The propensity scores are estimated using a
#' random forest. The `formula` supplied to `matchit()` is passed
#' directly to \pkgfun{randomForest}{randomForest}, and
#' \pkgfun{randomForest}{predict.randomForest} is used to compute the propensity
#' scores. The `link` argument is ignored, and predicted probabilities are
#' always returned as the distance measure. Note that because there is a random component, results will vary across runs unless a [seed][set.seed] is
#' set. }
#' \item{`"nnet"`}{ The
#' propensity scores are estimated using a single-hidden-layer neural network.
#' The `formula` supplied to `matchit()` is passed directly to
#' \pkgfun{nnet}{nnet}, and [fitted()] is used to compute the propensity scores.
#' The `link` argument is ignored, and predicted probabilities are always
#' returned as the distance measure. An argument to `size` must be
#' supplied to `distance.options` when using `method = "nnet"`. }
#' \item{`"cbps"`}{ The propensity scores are estimated using the
#' covariate balancing propensity score (CBPS) algorithm, which is a form of
#' logistic regression where balance constraints are incorporated to a
#' generalized method of moments estimation of of the model coefficients. The
#' `formula` supplied to `matchit()` is passed directly to
#' \pkgfun{CBPS}{CBPS}, and [fitted()] is used to compute the propensity
#' scores. The `link` argument can be specified as `"linear"` to use
#' the linear predictor instead of the predicted probabilities. No other links
#' are allowed. The `estimand` argument supplied to `matchit()` will
#' be used to select the appropriate estimand for use in defining the balance
#' constraints, so no argument needs to be supplied to `ATT` in
#' `CBPS`. }
#' \item{`"bart"`}{ The propensity scores are estimated
#' using Bayesian additive regression trees (BART). The `formula` supplied
#' to `matchit()` is passed directly to \pkgfun{dbarts}{bart2},
#' and \pkgfun{dbarts}{fitted.bart} is used to compute the propensity
#' scores. The `link` argument can be specified as `"linear"` to use
#' the linear predictor instead of the predicted probabilities. When
#' `s.weights` is supplied to `matchit()`, it will not be passed to
#' `bart2` because the `weights` argument in `bart2` does not
#' correspond to sampling weights. Note that because there is a random component to choosing the tuning
#' parameter, results will vary across runs unless the `seed` argument is supplied to `distance.options`. Note that setting a seed using [set.seed()] is not sufficient to guarantee reproducibility unless single-threading is used. See \pkgfun{dbarts}{bart2} for details.}
#' }
#'
#' ## Methods for computing distances from covariates
#'
#' The following methods involve computing a distance matrix from the covariates
#' themselves without estimating a propensity score. Calipers on the distance
#' measure and common support restrictions cannot be used, and the `distance`
#' component of the output object will be empty because no propensity scores are
#' estimated. The `link` and `distance.options` arguments are ignored with these
#' methods. See the individual matching methods pages for whether these
#' distances are allowed and how they are used. Each of these distance measures
#' can also be calculated outside `matchit()` using its [corresponding
#' function][euclidean_dist].
#'
#' \describe{
#' \item{`"euclidean"`}{ The Euclidean distance is the raw
#' distance between units, computed as \deqn{d_{ij} = \sqrt{(x_i - x_j)(x_i - x_j)'}} It is sensitive to the scale of the covariates, so covariates with
#' larger scales will take higher priority. }
#' \item{`"scaled_euclidean"`}{
#'  The scaled Euclidean distance is the
#' Euclidean distance computed on the scaled (i.e., standardized) covariates.
#' This ensures the covariates are on the same scale. The covariates are
#' standardized using the pooled within-group standard deviations, computed by
#' treatment group-mean centering each covariate before computing the standard
#' deviation in the full sample.
#'  }
#' \item{`"mahalanobis"`}{ The Mahalanobis distance is computed as \deqn{d_{ij} = \sqrt{(x_i - x_j)\Sigma^{-1}(x_i - x_j)'}} where \eqn{\Sigma} is the pooled within-group
#' covariance matrix of the covariates, computed by treatment group-mean
#' centering each covariate before computing the covariance in the full sample.
#' This ensures the variables are on the same scale and accounts for the
#' correlation between covariates. }
#' \item{`"robust_mahalanobis"`}{ The
#' robust rank-based Mahalanobis distance is the Mahalanobis distance computed
#' on the ranks of the covariates with an adjustment for ties. It is described
#' in Rosenbaum (2010, ch. 8) as an alternative to the Mahalanobis distance
#' that handles outliers and rare categories better than the standard
#' Mahalanobis distance but is not affinely invariant. }
#' }
#'
#' To perform Mahalanobis distance matching *and* estimate propensity scores to
#' be used for a purpose other than matching, the `mahvars` argument should be
#' used along with a different specification to `distance`. See the individual
#' matching method pages for details on how to use `mahvars`.
#'
#' ## Distances supplied as a numeric vector or matrix
#'
#' `distance` can also be supplied as a numeric vector whose values will be
#' taken to function like propensity scores; their pairwise difference will
#' define the distance between units. This might be useful for supplying
#' propensity scores computed outside `matchit()` or resupplying `matchit()`
#' with propensity scores estimated previously without having to recompute them.
#'
#' `distance` can also be supplied as a matrix whose values represent the
#' pairwise distances between units. The matrix should either be a square, with
#' a row and column for each unit (e.g., as the output of a call to
#' `as.matrix(`[`dist`]`(.))`), or have as many rows as there are treated units
#' and as many columns as there are control units (e.g., as the output of a call
#' to [mahalanobis_dist()] or \pkgfun{optmatch}{match_on}). Distance values of
#' `Inf` will disallow the corresponding units to be matched. When `distance` is
#' a supplied as a numeric vector or matrix, `link` and `distance.options` are
#' ignored.
#'
#' @note In versions of *MatchIt* prior to 4.0.0, `distance` was specified in a
#' slightly different way. When specifying arguments using the old syntax, they
#' will automatically be converted to the corresponding method in the new syntax
#' but a warning will be thrown. `distance = "logit"`, the old default, will
#' still work in the new syntax, though `distance = "glm", link = "logit"` is
#' preferred (note that these are the default settings and don't need to be made
#' explicit).
#'
#' @examples
#' data("lalonde")
#'
#' # Linearized probit regression PS:
#' m.out1 <- matchit(treat ~ age + educ + race + married +
#'                     nodegree + re74 + re75, data = lalonde,
#'                   distance = "glm", link = "linear.probit")
#' @examplesIf requireNamespace("mgcv", quietly = TRUE)
#' # GAM logistic PS with smoothing splines (s()):
#' m.out2 <- matchit(treat ~ s(age) + s(educ) + race + married +
#'                     nodegree + re74 + re75, data = lalonde,
#'                   distance = "gam")
#' summary(m.out2$model)
#' @examplesIf requireNamespace("CBPS", quietly = TRUE)
#' # CBPS for ATC matching w/replacement, using the just-
#' # identified version of CBPS (setting method = "exact"):
#' m.out3 <- matchit(treat ~ age + educ + race + married +
#'                     nodegree + re74 + re75, data = lalonde,
#'                   distance = "cbps", estimand = "ATC",
#'                   distance.options = list(method = "exact"),
#'                   replace = TRUE)
#' @examples
#' # Mahalanobis distance matching - no PS estimated
#' m.out4 <- matchit(treat ~ age + educ + race + married +
#'                     nodegree + re74 + re75, data = lalonde,
#'                   distance = "mahalanobis")
#'
#' m.out4$distance #NULL
#'
#' # Mahalanobis distance matching with PS estimated
#' # for use in a caliper; matching done on mahvars
#' m.out5 <- matchit(treat ~ age + educ + race + married +
#'                     nodegree + re74 + re75, data = lalonde,
#'                   distance = "glm", caliper = .1,
#'                   mahvars = ~ age + educ + race + married +
#'                                 nodegree + re74 + re75)
#'
#' summary(m.out5)
#'
#' # User-supplied propensity scores
#' p.score <- fitted(glm(treat ~ age + educ + race + married +
#'                         nodegree + re74 + re75, data = lalonde,
#'                       family = binomial))
#'
#' m.out6 <- matchit(treat ~ age + educ + race + married +
#'                     nodegree + re74 + re75, data = lalonde,
#'                   distance = p.score)
#'
#' # User-supplied distance matrix using optmatch::match_on()
#' @examplesIf requireNamespace("optmatch", quietly = TRUE)
#' dist_mat <- optmatch::match_on(
#'               treat ~ age + educ + race + nodegree +
#'                 married + re74 + re75, data = lalonde,
#'               method = "rank_mahalanobis")
#'
#' m.out7 <- matchit(treat ~ age + educ + race + nodegree +
#'                     married + re74 + re75, data = lalonde,
#'                   distance = dist_mat)

NULL

#distance2glm-----------------
distance2glm <- function(formula, data = NULL, link = "logit", ...) {

  linear <- is_not_null(link) && startsWith(as.character(link), "linear")
  if (linear) link <- sub("linear.", "", as.character(link), fixed = TRUE)

  args <- unique(c(names(formals(glm)), names(formals(glm.control))))
  A <- ...mget(args)
  A[lengths(A) == 0L] <- NULL

  family <- {
    if (is_null(A[["weights"]])) binomial(link = link)
    else quasibinomial(link = link)
  }

  A[["data"]] <- data
  A[["formula"]] <- formula
  A[["family"]] <- family

  res <- do.call("glm", A)

  pred <- predict(res, type = if (linear) "link" else "response")

  list(model = res, distance = pred)
}

#distance2gam-----------------
distance2gam <- function(formula, data = NULL, link = "logit", ...) {
  rlang::check_installed("mgcv")

  linear <- is_not_null(link) && startsWith(as.character(link), "linear")
  if (linear) link <- sub("linear.", "", as.character(link), fixed = TRUE)

  A <- list(...)
  weights <- A$weights
  A$weights <- NULL

  res <- do.call(mgcv::gam, c(list(formula, data, family = quasibinomial(link),
                                   weights = weights), A),
                 quote = TRUE)

  pred <- predict(res, type = if (linear) "link" else "response")

  list(model = res, distance = as.numeric(pred))
}

#distance2rpart-----------------
distance2rpart <- function(formula, data = NULL, link = NULL, ...) {
  rlang::check_installed("rpart")

  args <- unique(c(names(formals(rpart::rpart)), names(formals(rpart::rpart.control))))
  A <- ...mget(args)
  A[lengths(A) == 0L] <- NULL

  A$formula <- formula
  A$data <- data
  A$method <- "class"

  res <- do.call(rpart::rpart, A)
  list(model = res, distance = predict(res, type = "prob")[,"1"])
}

#distance2nnet-----------------
distance2nnet <- function(formula, data = NULL, link = NULL, ...) {
  rlang::check_installed("nnet")

  A <- list(...)
  weights <- A$weights
  A$weights <- NULL

  res <- do.call(nnet::nnet, c(list(formula, data = data, weights = weights, entropy = TRUE), A), quote = TRUE)
  list(model = res, distance = drop(fitted(res)))
}

#distance2cbps-----------------
distance2cbps <- function(formula, data = NULL, link = NULL, ...) {
  rlang::check_installed("CBPS")

  linear <- is_not_null(link) && startsWith(as.character(link), "linear")

  A <- list(...)

  A[["standardized"]] <- FALSE

  if (is_null(A[["ATT"]])) {
    if (is_null(A[["estimand"]])) {
      A[["ATT"]] <- 1
    }
    else {
      estimand <- toupper(A[["estimand"]])
      estimand <- match_arg(estimand, c("ATT", "ATC", "ATE"))
      A[["ATT"]] <- switch(estimand, "ATT" = 1, "ATC" = 2, 0)
    }
  }

  if (is_null(A[["method"]])) {
    A[["method"]] <- if (isFALSE(A[["over"]])) "exact" else "over"
  }

  A[c("estimand", "over")] <- NULL

  if (is_not_null(A[["weights"]])) {
    A[["sample.weights"]] <- A[["weights"]]
    A[["weights"]] <- NULL
  }

  A[["formula"]] <- formula
  A[["data"]] <- data

  capture.output({ #Keeps from printing message about treatment
    res <- do.call(CBPS::CBPS, A, quote = TRUE)
  })

  pred <- fitted(res)
  if (linear) pred <- qlogis(pred)

  list(model = res, distance = pred)
}

#distance2bart----------------
distance2bart <- function(formula, data = NULL, link = NULL, ...) {
  rlang::check_installed("dbarts")

  linear <- is_not_null(link) && startsWith(as.character(link), "linear")

  args <- unique(c(names(formals(dbarts::bart2)), names(formals(dbarts::dbartsControl))))
  A <- ...mget(args)
  A[lengths(A) == 0L] <- NULL

  A$formula <- formula
  A$data <- data

  res <- do.call(dbarts::bart2, A)

  pred <- fitted(res, type = if (linear) "link" else "response")

  list(model = res, distance = pred)
}

# distance2bart <- function(formula, data, link = NULL, ...) {
#   rlang::check_installed("BART")
#
#   if (is_not_null(link) && startsWith(as.character(link), "linear")) {
#     linear <- TRUE
#     link <- sub("linear.", "", as.character(link), fixed = TRUE)
#   }
#   else linear <- FALSE
#
#   #Keep link probit because default in matchit is logit but probit is much faster with BART
#   link <- "probit"
#
#   # if (is_null(link)) link <- "probit"
#   # else if (!link %in% c("probit", "logit")) {
#   #   stop("'link' must be \"probit\" or \"logit\" with distance = \"bart\".", call. = FALSE)
#   # }
#
#   data <- model.frame(formula, data)
#
#   treat <- binarize(data[[1]])
#   X <- data[-1]
#
#   chars <- vapply(X, is.character, logical(1L))
#   X[chars] <- lapply(X[chars], factor)
#
#   A <- list(...)
#
#   if (is_not_null(A[["mc.cores"]]) && A[["mc.cores"]][1] > 1) fun <- BART::mc.gbart
#   else fun <- BART::gbart
#
#   res <- do.call(fun, c(list(X,
#                              y.train = treat,
#                              type = switch(link, "logit" = "lbart", "pbart")),
#                         A[intersect(names(A), setdiff(names(formals(fun)),
#                                                       c("x.train", "y.train", "x.test", "type", "ntype")))]))
#
#   pred <- res$prob.train.mean
#   if (linear) pred <- switch(link, logit = qlogis, probit = qnorm)(pred)
#
#   return(list(model = res, distance = pred))
# }

#distance2randomforest-----------------
distance2randomforest <- function(formula, data = NULL, link = NULL, ...) {
  rlang::check_installed("randomForest")
  newdata <- get_all_vars(formula, data)
  treatvar <- as.character(formula[[2L]])
  newdata[[treatvar]] <- factor(newdata[[treatvar]], levels = c("0", "1"))
  res <- randomForest::randomForest(formula, data = newdata, ...)

  list(model = res, distance = predict(res, type = "prob")[,"1"])
}

#distance2glmnet--------------
distance2elasticnet <- function(formula, data = NULL, link = NULL, ...) {
  rlang::check_installed("glmnet")

  linear <- is_not_null(link) && startsWith(as.character(link), "linear")
  if (linear) link <- sub("linear.", "", as.character(link), fixed = TRUE)

  s <- ...get("s")
  if (is_null(s)) {
    s <- "lambda.1se"
  }

  args <- unique(c(names(formals(glmnet::glmnet)), names(formals(glmnet::cv.glmnet))))
  A <- ...mget(args)
  A[lengths(A) == 0L] <- NULL

  if (is_null(link)) link <- "logit"

  A$family <- switch(link,
                     "logit" = "binomial",
                     "log" = "poisson",
                     binomial(link = link))

  if (is_null(A[["alpha"]])) {
    A[["alpha"]] <- .5
  }

  mf <- model.frame(formula, data = data)

  A$y <- model.response(mf)
  A$x <- model.matrix(update(formula, . ~ . + 1), mf)[,-1,drop = FALSE]

  res <- do.call(glmnet::cv.glmnet, A)

  pred <- drop(predict(res, newx = A$x, s = s,
                       type = if (linear) "link" else "response"))

  list(model = res, distance = pred)
}
distance2lasso <- function(formula, data = NULL, link = NULL, ...) {
  if ("alpha" %in% ...names()) {
    args <- unique(c("s", names(formals(glmnet::glmnet)), names(formals(glmnet::cv.glmnet))))
    A <- ...mget(args)
    A[lengths(A) == 0L] <- NULL

    A$alpha <- 1
    do.call("distance2elasticnet", c(list(formula, data = data, link = link), A),
            quote = TRUE)
  }
  else {
    distance2elasticnet(formula = formula, data = data, link = link, alpha = 1, ...)
  }
}
distance2ridge <- function(formula, data = NULL, link = NULL, ...) {
  if ("alpha" %in% ...names()) {
    args <- unique(c("s", names(formals(glmnet::glmnet)), names(formals(glmnet::cv.glmnet))))
    A <- ...mget(args)
    A[lengths(A) == 0L] <- NULL

    A$alpha <- 0
    do.call("distance2elasticnet", c(list(formula, data = data, link = link), A),
            quote = TRUE)
  }
  else {
    distance2elasticnet(formula = formula, data = data, link = link, alpha = 0, ...)
  }
}

#distance2gbm--------------
distance2gbm <- function(formula, data = NULL, link = NULL, ...) {
  rlang::check_installed("gbm")

  linear <- is_not_null(link) && startsWith(as.character(link), "linear")

  method <- ...get("method")

  args <- names(formals(gbm::gbm))
  A <- ...mget(args)
  A[lengths(A) == 0L] <- NULL

  A$formula <- formula
  A$data <- data
  A$distribution <- "bernoulli"

  if (is_null(A[["n.trees"]])) A[["n.trees"]] <- 1e4
  if (is_null(A[["interaction.depth"]])) A[["interaction.depth"]] <- 3
  if (is_null(A[["shrinkage"]])) A[["shrinkage"]] <- .01
  if (is_null(A[["bag.fraction"]])) A[["bag.fraction"]] <- 1
  if (is_null(A[["cv.folds"]])) A[["cv.folds"]] <- 5
  if (is_null(A[["keep.data"]])) A[["keep.data"]] <- FALSE

  if (A[["cv.folds"]] <= 1 && A[["bag.fraction"]] == 1) {
    .err('either `bag.fraction` must be less than 1 or `cv.folds` must be greater than 1 when using `distance = "gbm"`')
  }

  if (is_null(method)) {
    if (A[["bag.fraction"]] < 1) method <- "OOB"
    else method <- "cv"
  }
  else if (!tolower(method) %in% c("oob", "cv")) {
    .err('`distance.options$method` should be one of "OOB" or "cv"')
  }

  res <- do.call(gbm::gbm, A)

  best.tree <- gbm::gbm.perf(res, plot.it = FALSE, method = method)

  pred <- drop(predict(res, newdata = data, n.trees = best.tree,
                       type = if (linear) "link" else "response"))

  list(model = res, distance = pred)
}
