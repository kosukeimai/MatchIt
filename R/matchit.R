#' Matching for Causal Inference
#'
#' @description
#' `matchit()` is the main function of *MatchIt* and performs
#' pairing, subset selection, and subclassification with the aim of creating
#' treatment and control groups balanced on included covariates. *MatchIt*
#' implements the suggestions of Ho, Imai, King, and Stuart (2007) for
#' improving parametric statistical models by preprocessing data with
#' nonparametric matching methods.
#'
#' This page documents the overall use of `matchit()`, but for specifics
#' of how `matchit()` works with individual matching methods, see the
#' individual pages linked in the Details section below.
#'
#' @param formula a two-sided [`formula`] object containing the treatment and
#' covariates to be used in creating the distance measure used in the matching.
#' This formula will be supplied to the functions that estimate the distance
#' measure. The formula should be specified as `A ~ X1 + X2 + ...` where
#' `A` represents the treatment variable and `X1` and `X2` are
#' covariates.
#' @param data a data frame containing the variables named in `formula`
#' and possible other arguments. If not found in `data`, the variables
#' will be sought in the environment.
#' @param method the matching method to be used. The allowed methods are
#' [`"nearest"`][method_nearest] for nearest neighbor matching (on
#' the propensity score by default), [`"optimal"`][method_optimal]
#' for optimal pair matching, [`"full"`][method_full] for optimal
#' full matching, [`"quick"`][method_quick] for generalized (quick)
#' full matching, [`"genetic"`][method_genetic] for genetic
#' matching, [`"cem"`][method_cem] for coarsened exact matching,
#' [`"exact"`][method_exact] for exact matching,
#' [`"cardinality"`][method_cardinality] for cardinality and
#' profile matching, and [`"subclass"`][method_subclass] for
#' subclassification. When set to `NULL`, no matching will occur, but
#' propensity score estimation and common support restrictions will still occur
#' if requested. See the linked pages for each method for more details on what
#' these methods do, how the arguments below are used by each on, and what
#' additional arguments are allowed.
#' @param distance the distance measure to be used. Can be either the name of a
#' method of estimating propensity scores (e.g., `"glm"`), the name of a
#' method of computing a distance matrix from the covariates (e.g.,
#' `"mahalanobis"`), a vector of already-computed distance measures, or a
#' matrix of pairwise distances. See [`distance`] for allowable
#' options. The default is `"glm"` for propensity scores estimated with
#' logistic regression using [glm()]. Ignored for some methods; see individual
#' methods pages for information on whether and how the distance measure is
#' used.
#' @param link when `distance` is specified as a string, an additional
#' argument controlling the link function used in estimating the distance
#' measure. Allowable options depend on the specific `distance` value
#' specified. See [`distance`] for allowable options with each
#' option. The default is `"logit"`, which, along with `distance = "glm"`, identifies the default measure as logistic regression propensity scores.
#' @param distance.options a named list containing additional arguments
#' supplied to the function that estimates the distance measure as determined
#' by the argument to `distance`. See [`distance`] for an
#' example of its use.
#' @param estimand a string containing the name of the target estimand desired.
#' Can be one of `"ATT"`, `"ATC"`, or `"ATE"`. Default is `"ATT"`. See Details and the individual methods
#' pages for information on how this argument is used.
#' @param exact for methods that allow it, for which variables exact matching
#' should take place. Can be specified as a string containing the names of
#' variables in `data` to be used or a one-sided formula with the desired
#' variables on the right-hand side (e.g., `~ X3 + X4`). See the
#' individual methods pages for information on whether and how this argument is
#' used.
#' @param mahvars for methods that allow it, on which variables Mahalanobis
#' distance matching should take place when `distance` corresponds to
#' propensity scores. Usually used to perform Mahalanobis distance matching
#' within propensity score calipers, where the propensity scores are computed
#' using `formula` and `distance`. Can be specified as a string
#' containing the names of variables in `data` to be used or a one-sided
#' formula with the desired variables on the right-hand side (e.g., `~ X3 + X4`). See the individual methods pages for information on whether and how this argument is used.
#' @param antiexact for methods that allow it, for which variables anti-exact
#' matching should take place. Anti-exact matching ensures paired individuals
#' do not have the same value of the anti-exact matching variable(s). Can be
#' specified as a string containing the names of variables in `data` to be
#' used or a one-sided formula with the desired variables on the right-hand
#' side (e.g., `~ X3 + X4`). See the individual methods pages for
#' information on whether and how this argument is used.
#' @param discard a string containing a method for discarding units outside a
#' region of common support. When a propensity score is estimated or supplied
#' to `distance` as a vector, the options are `"none"`,
#' `"treated"`, `"control"`, or `"both"`. For `"none"`, no
#' units are discarded for common support. Otherwise, units whose propensity
#' scores fall outside the corresponding region are discarded. Can also be a
#' `logical` vector where `TRUE` indicates the unit is to be
#' discarded. Default is `"none"` for no common support restriction. See
#' Details.
#' @param reestimate if `discard` is not `"none"` and propensity
#' scores are estimated, whether to re-estimate the propensity scores in the
#' remaining sample. Default is `FALSE` to use the propensity scores
#' estimated in the original sample.
#' @param s.weights an optional numeric vector of sampling weights to be
#' incorporated into propensity score models and balance statistics. Can also
#' be specified as a string containing the name of variable in `data` to
#' be used or a one-sided formula with the variable on the right-hand side
#' (e.g., `~ SW`). Not all propensity score models accept sampling
#' weights; see [`distance`] for information on which do and do not,
#' and see `vignette("sampling-weights")` for details on how to use
#' sampling weights in a matching analysis.
#' @param replace for methods that allow it, whether matching should be done
#' with replacement (`TRUE`), where control units are allowed to be
#' matched to several treated units, or without replacement (`FALSE`),
#' where control units can only be matched to one treated unit each. See the
#' individual methods pages for information on whether and how this argument is
#' used. Default is `FALSE` for matching without replacement.
#' @param m.order for methods that allow it, the order that the matching takes
#' place. Allowable options depend on the matching method. The default of
#' `NULL` corresponds to `"largest"` when a propensity score is
#' estimated or supplied as a vector and `"data"` otherwise.
#' @param caliper for methods that allow it, the width(s) of the caliper(s) to
#' use in matching. Should be a numeric vector with each value named according
#' to the variable to which the caliper applies. To apply to the distance
#' measure, the value should be unnamed. See the individual methods pages for
#' information on whether and how this argument is used. Positive values require the distance between paired units to be no larger than the supplied caliper; negative values require the distance between paired units to be larger than the absolute value value of the supplied caliper. The default is `NULL` for no caliper.
#' @param std.caliper `logical`; when a caliper is specified, whether the
#' the caliper is in standard deviation units (`TRUE`) or raw units
#' (`FALSE`). Can either be of length 1, applying to all calipers, or of
#' length equal to the length of `caliper`. Default is `TRUE`.
#' @param ratio for methods that allow it, how many control units should be
#' matched to each treated unit in k:1 matching. Should be a single integer
#' value. See the individual methods pages for information on whether and how
#' this argument is used. The default is 1 for 1:1 matching.
#' @param verbose `logical`; whether information about the matching
#' process should be printed to the console. What is printed depends on the
#' matching method. Default is `FALSE` for no printing other than
#' warnings.
#' @param include.obj `logical`; whether to include any objects created in
#' the matching process in the output, i.e., by the functions from other
#' packages `matchit()` calls. What is included depends on the matching
#' method. Default is `FALSE`.
#' @param normalize `logical`; whether to rescale the nonzero weights in each treatment group to have an average of 1. Default is `TRUE`. See "How Matching Weights Are Computed" below for more details.
#' @param \dots additional arguments passed to the functions used in the
#' matching process. See the individual methods pages for information on what
#' additional arguments are allowed for each method.
#'
#' @details
#' Details for the various matching methods can be found at the following help
#' pages:
#' * [`method_nearest`] for nearest neighbor matching
#' * [`method_optimal`] for optimal pair matching
#' * [`method_full`] for optimal full matching
#' * [`method_quick`] for generalized (quick) full matching
#' * [`method_genetic`] for genetic matching
#' * [`method_cem`] for coarsened exact matching
#' * [`method_exact`] for exact matching
#' * [`method_cardinality`] for cardinality and profile matching
#' * [`method_subclass`] for subclassification
#'
#' The pages contain information on what the method does, which of the arguments above are
#' allowed with them and how they are interpreted, and what additional
#' arguments can be supplied to further tune the method. Note that the default
#' method with no arguments supplied other than `formula` and `data`
#' is 1:1 nearest neighbor matching without replacement on a propensity score
#' estimated using a logistic regression of the treatment on the covariates.
#' This is not the same default offered by other matching programs, such as
#' those in *Matching*, `teffects` in Stata, or `PROC PSMATCH`
#' in SAS, so care should be taken if trying to replicate the results of those
#' programs.
#'
#' When `method = NULL`, no matching will occur, but any propensity score
#' estimation and common support restriction will. This can be a simple way to
#' estimate the propensity score for use in future matching specifications
#' without having to re-estimate it each time. The `matchit()` output with
#' no matching can be supplied to `summary()` to examine balance prior to
#' matching on any of the included covariates and on the propensity score if
#' specified. All arguments other than `distance`, `discard`, and
#' `reestimate` will be ignored.
#'
#' See [`distance`] for details on the several ways to
#' specify the `distance`, `link`, and `distance.options`
#' arguments to estimate propensity scores and create distance measures.
#'
#' When the treatment variable is not a `0/1` variable, it will be coerced
#' to one and returned as such in the `matchit()` output (see section
#' Value, below). The following rules are used: 1) if `0` is one of the
#' values, it will be considered the control and the other value the treated;
#' 2) otherwise, if the variable is a factor, `levels(treat)[1]` will be
#' considered control and the other value the treated; 3) otherwise,
#' `sort(unique(treat))[1]` will be considered control and the other value
#' the treated. It is safest to ensure the treatment variable is a `0/1`
#' variable.
#'
#' The `discard` option implements a common support restriction. It can
#' only be used when a distance measure is an estimated propensity score or supplied as a vector and is ignored for some matching
#' methods. When specified as `"treated"`, treated units whose distance
#' measure is outside the range of distance measures of the control units will
#' be discarded. When specified as `"control"`, control units whose
#' distance measure is outside the range of distance measures of the treated
#' units will be discarded. When specified as `"both"`, treated and
#' control units whose distance measure is outside the intersection of the
#' range of distance measures of the treated units and the range of distance
#' measures of the control units will be discarded. When `reestimate = TRUE` and `distance` corresponds to a propensity score-estimating
#' function, the propensity scores are re-estimated in the remaining units
#' prior to being used for matching or calipers.
#'
#' Caution should be used when interpreting effects estimated with various
#' values of `estimand`. Setting `estimand = "ATT"` doesn't
#' necessarily mean the average treatment effect in the treated is being
#' estimated; it just means that for matching methods, treated units will be
#' untouched and given weights of 1 and control units will be matched to them
#' (and the opposite for `estimand = "ATC"`). If a caliper is supplied or
#' treated units are removed for common support or some other reason (e.g.,
#' lacking matches when using exact matching), the actual estimand targeted is
#' not the ATT but the treatment effect in the matched sample. The argument to
#' `estimand` simply triggers which units are matched to which, and for
#' stratification-based methods (exact matching, CEM, full matching, and
#' subclassification), determines the formula used to compute the
#' stratification weights.
#'
#' ## How Matching Weights Are Computed
#'
#' Matching weights are computed in one of two ways depending on whether matching was done with replacement
#' or not.
#'
#' ### Matching without replacement and subclassification
#'
#' For matching *without* replacement (except for cardinality matching), including subclassification, each
#' unit is assigned to a subclass, which represents the pair they are a part of
#' (in the case of k:1 matching) or the stratum they belong to (in the case of
#' exact matching, coarsened exact matching, full matching, or
#' subclassification). The formula for computing the weights depends on the
#' argument supplied to `estimand`. A new "stratum propensity score"
#' (\eqn{p^s_i}) is computed for each unit \eqn{i} as \eqn{p^s_i = \frac{1}{n_s}\sum_{j: s_j =s_i}{I(A_j=1)}} where \eqn{n_s} is the size of subclass \eqn{s} and \eqn{I(A_j=1)} is 1 if unit \eqn{j} is treated and 0 otherwise. That is, the stratum propensity score for stratum \eqn{s} is the proportion of units in stratum \eqn{s} that are
#' in the treated group, and all units in stratum \eqn{s} are assigned that
#' stratum propensity score. This is distinct from the propensity score used for matching, if any. Weights are then computed using the standard formulas for
#' inverse probability weights with the stratum propensity score inserted:
#' * for the ATT, weights are 1 for the treated
#' units and \eqn{\frac{p^s}{1-p^s}} for the control units
#' * for the ATC, weights are
#' \eqn{\frac{1-p^s}{p^s}} for the treated units and 1 for the control units
#' * for the ATE, weights are \eqn{\frac{1}{p^s}} for the treated units and \eqn{\frac{1}{1-p^s}} for the
#' control units.
#'
#' For cardinality matching, all matched units receive a weight
#' of 1.
#'
#' ### Matching with replacement
#'
#' For matching *with* replacement, units are not assigned to unique strata. For
#' the ATT, each treated unit gets a weight of 1. Each control unit is weighted
#' as the sum of the inverse of the number of control units matched to the same
#' treated unit across its matches. For example, if a control unit was matched
#' to a treated unit that had two other control units matched to it, and that
#' same control was matched to a treated unit that had one other control unit
#' matched to it, the control unit in question would get a weight of \eqn{1/3 + 1/2 = 5/6}. For the ATC, the same is true with the treated and control labels
#' switched. The weights are computed using the `match.matrix` component
#' of the `matchit()` output object.
#'
#' ### Normalized weights
#'
#' When `normalize = TRUE` (the default), in each treatment group, weights are divided by the mean of the nonzero
#' weights in that treatment group to make the weights sum to the number of
#' units in that treatment group (i.e., to have an average of 1).
#'
#' ### Sampling weights
#'
#' If sampling weights are included through the
#' `s.weights` argument, they will be included in the `matchit()`
#' output object but not incorporated into the matching weights.
#' [match_data()], which extracts the matched set from a `matchit` object,
#' combines the matching weights and sampling weights.
#'
#' @return When `method` is something other than `"subclass"`, a
#' `matchit` object with the following components:
#'
#' \item{match.matrix}{a matrix containing the matches. The row names correspond
#' to the treated units and the values in each row are the names (or indices)
#' of the control units matched to each treated unit. When treated units are
#' matched to different numbers of control units (e.g., with variable ratio matching or
#' matching with a caliper), empty spaces will be filled with `NA`. Not
#' included when `method` is `"full"`, `"cem"` (unless `k2k = TRUE`), `"exact"`, `"quick"`, or `"cardinality"` (unless `mahvars` is supplied and `ratio` is an integer).}
#' \item{subclass}{a factor
#' containing matching pair/stratum membership for each unit. Unmatched units
#' will have a value of `NA`. Not included when `replace = TRUE` or when `method = "cardinality"` unless `mahvars` is supplied and `ratio` is an integer.}
#' \item{weights}{a numeric vector of estimated matching weights. Unmatched and
#' discarded units will have a weight of zero.}
#' \item{model}{the fit object of
#' the model used to estimate propensity scores when `distance` is
#' specified as a method of estimating propensity scores. When
#' `reestimate = TRUE`, this is the model estimated after discarding
#' units.}
#' \item{X}{a data frame of covariates mentioned in `formula`, `exact`, `mahvars`, `caliper`, and `antiexact`.}
#' \item{call}{the `matchit()` call.}
#' \item{info}{information on the matching method and distance measures used.}
#' \item{estimand}{the argument supplied to `estimand`.}
#' \item{formula}{the `formula` supplied.}
#' \item{treat}{a vector of treatment status converted to zeros (0) and ones
#' (1) if not already in that format.}
#' \item{distance}{a vector of distance
#' values (i.e., propensity scores) when `distance` is supplied as a
#' method of estimating propensity scores or a numeric vector.}
#' \item{discarded}{a logical vector denoting whether each observation was
#' discarded (`TRUE`) or not (`FALSE`) by the argument to `discard`.}
#' \item{s.weights}{the vector of sampling weights supplied to the `s.weights` argument, if any.}
#' \item{exact}{a one-sided formula containing the variables, if any, supplied to `exact`.}
#' \item{mahvars}{a one-sided formula containing the variables, if any, supplied to `mahvars`.}
#' \item{obj}{when `include.obj = TRUE`, an object containing the intermediate results of the matching procedure. See
#' the individual methods pages for what this component will contain.}
#'
#' When `method = "subclass"`, a `matchit.subclass` object with the same
#' components as above except that `match.matrix` is excluded and one
#' additional component, `q.cut`, is included, containing a vector of the
#' distance measure cutpoints used to define the subclasses. See
#' [`method_subclass`] for details.
#'
#' @author Daniel Ho, Kosuke Imai, Gary King, and Elizabeth Stuart wrote the original package. Starting with version 4.0.0, Noah Greifer is the primary maintainer and developer.
#'
#' @seealso [summary.matchit()] for balance assessment after matching, [plot.matchit()] for plots of covariate balance and propensity score overlap after matching.
#'
#' * `vignette("MatchIt")` for an introduction to matching with *MatchIt*
#' * `vignette("matching-methods")` for descriptions of the variety of matching methods and options available
#' * `vignette("assessing-balance")` for information on assessing the quality of a matching specification
#' * `vignette("estimating-effects")` for instructions on how to estimate treatment effects after matching
#' * `vignette("sampling-weights")` for a guide to using *MatchIt* with sampling weights.
#'
#' @references
#' Ho, D. E., Imai, K., King, G., & Stuart, E. A. (2007). Matching
#' as Nonparametric Preprocessing for Reducing Model Dependence in Parametric
#' Causal Inference. *Political Analysis*, 15(3), 199–236. \doi{10.1093/pan/mpl013}
#'
#' Ho, D. E., Imai, K., King, G., & Stuart, E. A. (2011). MatchIt:
#' Nonparametric Preprocessing for Parametric Causal Inference. *Journal of Statistical Software*, 42(8). \doi{10.18637/jss.v042.i08}
#'
#' @examples
#' data("lalonde")
#'
#' # Default: 1:1 NN PS matching w/o replacement
#'
#' m.out1 <- matchit(treat ~ age + educ + race + nodegree +
#'                     married + re74 + re75,
#'                   data = lalonde)
#' m.out1
#' summary(m.out1)
#'
#' # 1:1 NN Mahalanobis distance matching w/ replacement and
#' # exact matching on married and race
#'
#' m.out2 <- matchit(treat ~ age + educ + race + nodegree +
#'                     married + re74 + re75,
#'                   data = lalonde,
#'                   distance = "mahalanobis",
#'                   replace = TRUE,
#'                   exact = ~ married + race)
#' m.out2
#' summary(m.out2, un = TRUE)
#'
#' # 2:1 NN Mahalanobis distance matching within caliper defined
#' # by a probit pregression PS
#'
#' m.out3 <- matchit(treat ~ age + educ + race + nodegree +
#'                     married + re74 + re75,
#'                   data = lalonde,
#'                   distance = "glm",
#'                   link = "probit",
#'                   mahvars = ~ age + educ + re74 + re75,
#'                   caliper = .1,
#'                   ratio = 2)
#' m.out3
#' summary(m.out3, un = TRUE)
#'
#' # Optimal full PS matching for the ATE within calipers on
#' # PS, age, and educ
#' @examplesIf requireNamespace("optmatch", quietly = TRUE)
#' m.out4 <- matchit(treat ~ age + educ + race + nodegree +
#'                     married + re74 + re75,
#'                   data = lalonde,
#'                   method = "full",
#'                   estimand = "ATE",
#'                   caliper = c(.1, age = 2, educ = 1),
#'                   std.caliper = c(TRUE, FALSE, FALSE))
#' m.out4
#' summary(m.out4, un = TRUE)
#' @examples
#' # Subclassification on a logistic PS with 10 subclasses after
#' # discarding controls outside common support of PS
#'
#' s.out1 <- matchit(treat ~ age + educ + race + nodegree +
#'                     married + re74 + re75,
#'                   data = lalonde,
#'                   method = "subclass",
#'                   distance = "glm",
#'                   discard = "control",
#'                   subclass = 10)
#' s.out1
#' summary(s.out1, un = TRUE)

#' @export
matchit <- function(formula,
                    data = NULL,
                    method = "nearest",
                    distance = "glm",
                    link = "logit",
                    distance.options = list(),
                    estimand = "ATT",
                    exact = NULL,
                    mahvars = NULL,
                    antiexact = NULL,
                    discard = "none",
                    reestimate = FALSE,
                    s.weights = NULL,
                    replace = FALSE,
                    m.order = NULL,
                    caliper = NULL,
                    std.caliper = TRUE,
                    ratio = 1,
                    verbose = FALSE,
                    include.obj = FALSE,
                    normalize = TRUE,
                    ...) {

  #Checking input format
  #data input
  mcall <- match.call()

  ## Process method
  chk::chk_null_or(method, vld = chk::vld_string)
  if (is_null(method)) {
    fn2 <- "matchit2null"
  }
  else {
    method <- tolower(method)
    method <- match_arg(method, c("exact", "cem", "nearest", "optimal",
                                  "full", "genetic", "subclass", "cardinality",
                                  "quick"))
    fn2 <- paste0("matchit2", method)
  }

  #Process formula and data inputs
  if (!rlang::is_formula(formula, lhs = TRUE)) {
    .err("`formula` must be a formula relating treatment to covariates")
  }

  treat.form <- update(terms(formula, data = data), . ~ 0)
  treat.mf <- model.frame(treat.form, data = data, na.action = "na.pass")
  treat <- model.response(treat.mf)

  #Check and binarize treat
  treat <- check_treat(treat)
  if (is_null(treat)) {
    .err("the treatment cannot be `NULL`")
  }

  names(treat) <- rownames(treat.mf)

  n.obs <- length(treat)

  #Process inputs
  ignored.inputs <- check.inputs(mcall = mcall, method = method, distance = distance,
                                 link = link, distance.options = distance.options, exact = exact,
                                 mahvars = mahvars, antiexact = antiexact, caliper = caliper, discard = discard,
                                 reestimate = reestimate, s.weights = s.weights, replace = replace,
                                 ratio = ratio, m.order = m.order, estimand = estimand)

  for (i in ignored.inputs) {
    assign(i, NULL)
  }

  #Process replace
  replace <- process.replace(replace, method, ...)

  #Process ratio
  ratio <- process.ratio(ratio, method, ...)

  #Process s.weights
  if (is_not_null(s.weights)) {
    if (is.character(s.weights)) {
      if (is_null(data) || !is.data.frame(data)) {
        .err("if `s.weights` is specified a string, a data frame containing the named variable must be supplied to `data`")
      }

      if (!all(hasName(data, s.weights))) {
        .err("the name supplied to `s.weights` must be a variable in `data`")
      }

      s.weights.form <- reformulate(s.weights)
      s.weights <- model.frame(s.weights.form, data, na.action = "na.pass")

      if (ncol(s.weights) != 1L) {
        .err("`s.weights` can only contain one named variable")
      }

      s.weights <- s.weights[[1L]]
    }
    else if (rlang::is_formula(s.weights)) {
      s.weights.form <- update(terms(s.weights, data = data), NULL ~ .)
      s.weights <- model.frame(s.weights.form, data, na.action = "na.pass")

      if (ncol(s.weights) != 1L) {
        .err("`s.weights` can only contain one named variable")
      }

      s.weights <- s.weights[[1L]]
    }
    else if (!is.numeric(s.weights)) {
      .err("`s.weights` must be supplied as a numeric vector, string, or one-sided formula")
    }

    chk::chk_not_any_na(s.weights)
    if (length(s.weights) != n.obs) {
      .err("`s.weights` must be the same length as the treatment vector")
    }

    names(s.weights) <- names(treat)
  }

  #Process distance function
  is.full.mahalanobis <- FALSE
  fn1 <- NULL
  if (is_null(method) || !method %in% c("exact", "cem", "cardinality")) {
    distance <- process.distance(distance, method, treat)

    if (is.numeric(distance)) {
      fn1 <- "distance2user"
    }
    else if (is.character(distance)) {
      if (distance %in% matchit_distances()) {
        fn1 <- "distance2mahalanobis"
        is.full.mahalanobis <- TRUE
        attr(is.full.mahalanobis, "transform") <- distance
      }
      else {
        fn1 <- paste0("distance2", distance)
      }
    }
  }

  #Process covs
  if (is_not_null(fn1) && fn1 == "distance2gam") {
    rlang::check_installed("mgcv")
    env <- environment(formula)
    covs.formula <- mgcv::interpret.gam(formula)$fake.formula
    environment(covs.formula) <- env
    covs.formula <- delete.response(terms(covs.formula, data = data))
  }
  else {
    covs.formula <- delete.response(terms(formula, data = data))
  }

  covs.formula <- update(covs.formula, ~ .)
  covs <- model.frame(covs.formula, data = data, na.action = "na.pass")
  k <- ncol(covs)
  for (i in seq_len(k)) {
    if (anyNA(covs[[i]]) || (is.numeric(covs[[i]]) && !all(is.finite(covs[[i]])))) {
      covariates.with.missingness <- names(covs)[i:k][vapply(i:k, function(j) anyNA(covs[[j]]) ||
                                                               (is.numeric(covs[[j]]) && !all(is.finite(covs[[j]]))),
                                                             logical(1L))]
      .err(paste0("Missing and non-finite values are not allowed in the covariates. Covariates with missingness or non-finite values:\n\t",
                  toString(covariates.with.missingness)), tidy = FALSE)
    }

    if (is.character(covs[[i]])) {
      covs[[i]] <- factor(covs[[i]])
    }
  }

  #Process exact, mahvars, and antiexact
  exactcovs <- process.variable.input(exact, data)
  exact <- attr(exactcovs, "terms")

  mahcovs <- process.variable.input(mahvars, data)
  mahvars <- attr(mahcovs, "terms")

  antiexactcovs <- process.variable.input(antiexact, data)
  antiexact <- attr(antiexactcovs, "terms")

  chk::chk_flag(verbose)
  chk::chk_flag(normalize)

  #Estimate distance, discard from common support, optionally re-estimate distance
  if (is_null(fn1) || is.full.mahalanobis) {
    #No distance measure
    dist.model <- distance <- link <- NULL
  }
  else if (fn1 == "distance2user") {
    dist.model <- link <- NULL
  }
  else {
    .cat_verbose("Estimating propensity scores...\n", verbose = verbose)

    if (is_not_null(s.weights)) {
      attr(s.weights, "in_ps") <- !distance %in% c("bart")
    }

    #Estimate distance
    if (is_null(distance.options)) {
      distance.options <- list(formula = formula,
                               data = data,
                               verbose = verbose,
                               estimand = estimand)
    }
    else {
      chk::chk_list(distance.options)

      if (is_null(distance.options$formula)) distance.options$formula <- formula
      if (is_null(distance.options$data)) distance.options$data <- data
      if (is_null(distance.options$verbose)) distance.options$verbose <- verbose
      if (is_null(distance.options$estimand)) distance.options$estimand <- estimand
    }

    if (is_null(distance.options$weights) && !fn1 %in% c("distance2bart")) {
      distance.options$weights <- s.weights
    }

    distance.options$link <- {
      if (is_not_null(attr(distance, "link"))) attr(distance, "link")
      else link
    }

    dist.out <- do.call(fn1, distance.options)

    dist.model <- dist.out$model
    distance <- dist.out$distance

    #Remove smoothing terms from gam formula
    if (inherits(dist.model, "gam")) {
      env <- environment(formula)
      formula <- mgcv::interpret.gam(formula)$fake.formula
      environment(formula) <- env
    }
  }

  #Process discard
  if (is_null(fn1) || is.full.mahalanobis || identical(fn1, "distance2user")) {
    discarded <- discard(treat, distance, discard)
  }
  else {
    discarded <- discard(treat, dist.out$distance, discard)

    #Optionally reestimate
    if (reestimate && any(discarded)) {
      for (i in seq_along(distance.options)) {
        if (length(distance.options[[i]]) == n.obs) {
          distance.options[[i]] <- distance.options[[i]][!discarded]
        }
        else if (length(dim(distance.options[[i]])) == 2L && nrow(distance.options[[i]]) == n.obs) {
          distance.options[[i]] <- distance.options[[i]][!discarded, , drop = FALSE]
        }
      }
      dist.out <- do.call(fn1, distance.options, quote = TRUE)
      dist.model <- dist.out$model
      distance[!discarded] <- dist.out$distance
    }
  }

  #Process caliper
  calcovs <- NULL

  if (is_not_null(caliper)) {
    caliper <- process.caliper(caliper, method, data, covs, mahcovs, distance, discarded, std.caliper)

    if (is_not_null(attr(caliper, "cal.formula"))) {
      calcovs <- model.frame(attr(caliper, "cal.formula"), data = data,
                             na.action = "na.pass")

      if (anyNA(calcovs)) {
        .err("missing values are not allowed in the covariates named in `caliper`")
      }

      attr(caliper, "cal.formula") <- NULL
    }
  }

  #Matching!
  match.out <- do.call(fn2, list(treat = treat, covs = covs, data = data, distance = distance,
                                 discarded = discarded, exact = exact, mahvars = mahvars,
                                 replace = replace, m.order = m.order, caliper = caliper,
                                 s.weights = s.weights, ratio = ratio, is.full.mahalanobis = is.full.mahalanobis,
                                 formula = formula, estimand = estimand, verbose = verbose,
                                 antiexact = antiexact, ...),
                       quote = TRUE)

  weights <- match.out[["weights"]]

  #Normalize weights
  if (normalize) {
    wi <- which(weights > 0)
    weights[wi] <- .make_sum_to_n(weights[wi], treat[wi])
  }

  info <- create_info(method, fn1, link, discard, replace, ratio,
                      mahalanobis = is.full.mahalanobis || is_not_null(mahvars),
                      transform = attr(is.full.mahalanobis, "transform"),
                      subclass = match.out$subclass,
                      antiexact = colnames(antiexactcovs),
                      distance_is_matrix = is_not_null(distance) && is.matrix(distance))

  #Create X output, removing duplicate variables
  X.list.nm <- c("covs", "exactcovs", "mahcovs", "calcovs", "antiexactcovs")
  X <- NULL
  for (i in X.list.nm) {
    X_tmp <- get0(i, inherits = FALSE)

    if (is_null(X_tmp)) {
      next
    }

    if (is_null(X)) {
      X <- X_tmp
    }
    else if (!all(hasName(X, names(X_tmp)))) {
      X <- cbind(X, X_tmp[!names(X_tmp) %in% names(X)])
    }
  }

  ## putting all the results together
  out <- list(
    match.matrix = match.out[["match.matrix"]],
    subclass = match.out[["subclass"]],
    weights = weights,
    X = X,
    call = mcall,
    info = info,
    estimand = estimand,
    formula = formula,
    treat = treat,
    distance = if (is_not_null(distance) && !is.matrix(distance)) setNames(distance, names(treat)),
    discarded = discarded,
    s.weights = s.weights,
    exact = exact,
    mahvars = mahvars,
    caliper = caliper,
    q.cut = match.out[["q.cut"]],
    model = dist.model,
    obj = if (include.obj) match.out[["obj"]]
  )

  out[lengths(out) == 0L] <- NULL

  class(out) <- class(match.out)

  out
}

#' @exportS3Method print matchit
print.matchit <- function(x, ...) {
  info <- x[["info"]]
  cal <- is_not_null(x[["caliper"]])
  dis <- c("both", "control", "treat")[pmatch(info$discard, c("both", "control", "treat"), 0L)]
  disl <- is_not_null(dis)
  nm <- is_null(x[["method"]])

  cat("A `matchit` object\n")

  cat(sprintf(" - method: %s\n", info_to_method(info)))

  if (is_not_null(info$distance) || info$mahalanobis) {
    cat(" - distance: ")
    if (info$mahalanobis) {
      if (is_null(info$transform)) #mahvars used
        cat("Mahalanobis")
      else {
        cat(capwords(gsub("_", " ", info$transform, fixed = TRUE)))
      }
    }

    if (is_not_null(info$distance) && !info$distance %in% matchit_distances()) {
      if (info$mahalanobis) cat(" [matching]\n             ")

      if (info$distance_is_matrix) cat("User-defined (matrix)")
      else if (info$distance != "user") cat("Propensity score")
      else if (is_not_null(attr(info$distance, "custom"))) cat(attr(info$distance, "custom"))
      else cat("User-defined")

      if (cal || disl) {
        cal.ps <- hasName(x[["caliper"]], "")
        cat(sprintf(" [%s]\n",
                    toString(c("matching", "subclassification", "caliper", "common support")[c(!nm && !info$mahalanobis && info$method != "subclass", !nm && info$method == "subclass", cal.ps, disl)])))
      }

      if (info$distance != "user") {
        cat(sprintf("\n             - estimated with %s\n",
                    info_to_distance(info)))
        if (is_not_null(x[["s.weights"]])) {
          cat(sprintf("             - sampling weights %s in estimation\n",
                      if (isTRUE(attr(x[["s.weights"]], "in_ps"))) "included" else "not included"))
        }
      }
    }
  }

  if (cal) {
    cat(sprintf(" - caliper: %s\n",
                toString(vapply(seq_along(x[["caliper"]]),
                                function(z) {
                                  sprintf("%s (%s)",
                                          if (nzchar(names(x[["caliper"]])[z])) names(x[["caliper"]])[z]
                                          else "<distance>",
                                          format(round(x[["caliper"]][z], 3L)))
                                }, character(1L)))))
  }

  if (disl) {
    cat(sprintf(" - common support: %s dropped\n",
                switch(dis,
                       "both" = "units from both groups",
                       "treat" = "treated units",
                       "control" = "control units")))
  }

  cat(sprintf(" - number of obs.: %s (original)%s\n",
              length(x[["treat"]]),
              if (all_equal_to(x[["weights"]], 1)) ""
              else sprintf(", %s (matched)", sum(x[["weights"]] != 0))))

  if (is_not_null(x[["s.weights"]])) {
    cat(" - sampling weights: present\n")
  }

  if (is_not_null(x[["estimand"]])) {
    cat(sprintf(" - target estimand: %s\n", x[["estimand"]]))
  }

  if (is_not_null(x[["X"]])) {
    cat(sprintf(" - covariates: %s\n",
                if (length(names(x[["X"]])) > 40L) "too many to name"
                else toString(names(x[["X"]]))))
  }

  invisible(x)
}

matchit2null <- function(discarded, ...) {

  res <- list(weights = as.numeric(!discarded))
  class(res) <- "matchit"

  res
}
