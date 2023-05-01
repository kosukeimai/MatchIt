#' Matching for Causal Inference
#'
#' @aliases matchit print.matchit
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
#' full matching, [`"genetic"`][method_genetic] for genetic
#' matching, [`"cem"`][method_cem] for coarsened exact matching,
#' [`"exact"`][method_exact] for exact matching,
#' [`"cardinality"`][method_cardinality] for cardinality and
#' template matching, and [`"subclass"`][method_subclass] for
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
#' option. The default is `"logit"`, which, along with `distance = "glm"`, identifies the default measure as logistic regression propensity
#' scores.
#' @param distance.options a named list containing additional arguments
#' supplied to the function that estimates the distance measure as determined
#' by the argument to `distance`. See [distance] for an
#' example of its use.
#' @param estimand a string containing the name of the target estimand desired.
#' Can be one of `"ATT"` or `"ATC"`. Some methods accept `"ATE"`
#' as well. Default is `"ATT"`. See Details and the individual methods
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
#' formula with the desired variables on the right-hand side (e.g., `~ X3 + X4`). See the individual methods pages for information on whether and how
#' this argument is used.
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
#' weights; see [distance] for information on which do and do not,
#' and see `vignette("sampling-weights")` for details on how to use
#' sampling weights in a matching analysis.
#' @param replace for methods that allow it, whether matching should be done
#' with replacement (`TRUE`), where control units are allowed to be
#' matched to several treated units, or without replacement (`FALSE`),
#' where control units can only be matched to one treated unit each. See the
#' individual methods pages for information on whether and how this argument is
#' used. Default is `FALSE` for matching without replacement.
#' @param m.order for methods that allow it, the order that the matching takes
#' place. Allowable options depend on the matching method but include
#' `"largest"`, where matching takes place in descending order of distance
#' measures; `"smallest"`, where matching takes place in ascending order
#' of distance measures; `"random"`, where matching takes place in a
#' random order; and `"data"` where matching takes place based on the
#' order of units in the data. When `m.order = "random"`, results may
#' differ across different runs of the same code unless a seed is set and
#' specified with [set.seed()]. See the individual methods pages for
#' information on whether and how this argument is used. The default of
#' `NULL` corresponds to `"largest"` when a propensity score is
#' estimated or supplied as a vector and `"data"` otherwise.
#' @param caliper for methods that allow it, the width(s) of the caliper(s) to
#' use in matching. Should be a numeric vector with each value named according
#' to the variable to which the caliper applies. To apply to the distance
#' measure, the value should be unnamed. See the individual methods pages for
#' information on whether and how this argument is used. The default is
#' `NULL` for no caliper.
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
#' @param \dots additional arguments passed to the functions used in the
#' matching process. See the individual methods pages for information on what
#' additional arguments are allowed for each method. Ignored for `print()`.
#' @param x a `matchit` object.
#'
#' @details
#' Details for the various matching methods can be found at the following help
#' pages:
#' * [`method_nearest`] for nearest neighbor matching
#' * [`method_optimal`] for optimal pair matching
#' * [`method_full`] for optimal full matching
#' * [`method_genetic`] for genetic matching
#' * [`method_cem`] for coarsened exact matching
#' * [`method_exact`] for exact matching
#' * [`method_cardinality`] for cardinality and template matching
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
#' See [distance] for details on the several ways to
#' specify the `distance`, `link`, and `distance.options`
#' arguments to estimate propensity scores and create distance measures.
#'
#' When the treatment variable is not a `0/1` variable, it will be coerced
#' to one and returned as such in the `matchit()` output (see section
#' Value, below). The following rules are used: 1) if `0` is one of the
#' values, it will be considered the control and the other value the treated;
#' 2) otherwise, if the variable is a factor, `levels(treat)[1]` will be
#' considered control and the other variable the treated; 3) otherwise,
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
#' For matching *without* replacement (except for cardinality matching), each
#' unit is assigned to a subclass, which represents the pair they are a part of
#' (in the case of k:1 matching) or the stratum they belong to (in the case of
#' exact matching, coarsened exact matching, full matching, or
#' subclassification). The formula for computing the weights depends on the
#' argument supplied to `estimand`. A new stratum "propensity score"
#' (`p`) is computed as the proportion of units in each stratum that are
#' in the treated group, and all units in that stratum are assigned that
#' propensity score. Weights are then computed using the standard formulas for
#' inverse probability weights: for the ATT, weights are 1 for the treated
#' units and `p/(1-p)` for the control units; for the ATC, weights are
#' `(1-p)/p` for the treated units and 1 for the control units; for the
#' ATE, weights are `1/p` for the treated units and `1/(1-p)` for the
#' control units. For cardinality matching, all matched units receive a weight
#' of 1.
#'
#' For matching *with* replacement, units are not assigned to unique strata. For
#' the ATT, each treated unit gets a weight of 1. Each control unit is weighted
#' as the sum of the inverse of the number of control units matched to the same
#' treated unit across its matches. For example, if a control unit was matched
#' to a treated unit that had two other control units matched to it, and that
#' same control was matched to a treated unit that had one other control unit
#' matched to it, the control unit in question would get a weight of 1/3 + 1/2
#' = 5/6. For the ATC, the same is true with the treated and control labels
#' switched. The weights are computed using the `match.matrix` component
#' of the `matchit()` output object.
#'
#' In each treatment group, weights are divided by the mean of the nonzero
#' weights in that treatment group to make the weights sum to the number of
#' units in that treatment group. If sampling weights are included through the
#' `s.weights` argument, they will be included in the `matchit()`
#' output object but not incorporated into the matching weights.
#' [match.data()], which extracts the matched set from a `matchit` object,
#' combines the matching weights and sampling weights.
#'
#' @return When `method` is something other than `"subclass"`, a
#' `matchit` object with the following components:
#'
#' \item{match.matrix}{a matrix containing the matches. The rownames correspond
#' to the treated units and the values in each row are the names (or indices)
#' of the control units matched to each treated unit. When treated units are
#' matched to different numbers of control units (e.g., with exact matching or
#' matching with a caliper), empty spaces will be filled with `NA`. Not
#' included when `method` is `"full"`, `"cem"` (unless `k2k
#' = TRUE`), `"exact"`, or `"cardinality"`.}
#' \item{subclass}{a factor
#' containing matching pair/stratum membership for each unit. Unmatched units
#' will have a value of `NA`. Not included when `replace = TRUE`.}
#' \item{weights}{a numeric vector of estimated matching weights. Unmatched and
#' discarded units will have a weight of zero.}
#' \item{model}{the fit object of
#' the model used to estimate propensity scores when `distance` is
#' specified and not `"mahalanobis"` or a numeric vector. When
#' `reestimate = TRUE`, this is the model estimated after discarding
#' units.}
#' \item{X}{a data frame of covariates mentioned in `formula`,
#' `exact`, `mahvars`, and `antiexact`.}
#' \item{call}{the `matchit()` call.}
#' \item{info}{information on the matching method and
#' distance measures used.}
#' \item{estimand}{the argument supplied to
#' `estimand`.}
#' \item{formula}{the `formula` supplied.}
#' \item{treat}{a vector of treatment status converted to zeros (0) and ones
#' (1) if not already in that format.}
#' \item{distance}{a vector of distance
#' values (i.e., propensity scores) when `distance` is supplied as a
#' method of estimating propensity scores or a numeric vector.}
#' \item{discarded}{a logical vector denoting whether each observation was
#' discarded (`TRUE`) or not (`FALSE`) by the argument to
#' `discard`.}
#' \item{s.weights}{the vector of sampling weights supplied to
#' the `s.weights` argument, if any.}
#' \item{exact}{a one-sided formula
#' containing the variables, if any, supplied to `exact`.}
#' \item{mahvars}{a one-sided formula containing the variables, if any,
#' supplied to `mahvars`.}
#' \item{obj}{when `include.obj = TRUE`, an
#' object containing the intermediate results of the matching procedure. See
#' the individual methods pages for what this component will contain.}
#'
#' When `method = "subclass"`, a `matchit.subclass` object with the same
#' components as above except that `match.matrix` is excluded and one
#' additional component, `q.cut`, is included, containing a vector of the
#' distance measure cutpoints used to define the subclasses. See
#' [`method_subclass`] for details.
#'
#' @author Daniel Ho (\email{dho@@law.stanford.edu}); Kosuke Imai
#' (\email{imai@@harvard.edu}); Gary King (\email{king@@harvard.edu});
#' Elizabeth Stuart (\email{estuart@@jhsph.edu})
#'
#' Version 4.0.0 update by Noah Greifer (\email{noah.greifer@@gmail.com})
#'
#' @seealso [summary.matchit()] for balance assessment after matching, [plot.matchit()] for plots of covariate balance and propensity score overlap after matching.
#'
#' `vignette("MatchIt")` for an introduction to matching with
#' *MatchIt*; `vignette("matching-methods")` for descriptions of the
#' variety of matching methods and options available;
#' `vignette("assessing-balance")` for information on assessing the
#' quality of a matching specification; `vignette("estimating-effects")`
#' for instructions on how to estimate treatment effects after matching; and
#' `vignette("sampling-weights")` for a guide to using *MatchIt* with
#' sampling weights.
#'
#' @references Ho, D. E., Imai, K., King, G., & Stuart, E. A. (2007). Matching
#' as Nonparametric Preprocessing for Reducing Model Dependence in Parametric
#' Causal Inference. *Political Analysis*, 15(3), 199–236. \doi{10.1093/pan/mpl013}
#'
#' Ho, D. E., Imai, K., King, G., & Stuart, E. A. (2011). MatchIt:
#' Nonparametric Preprocessing for Parametric Causal Inference. *Journal of
#' Statistical Software*, 42(8). \doi{10.18637/jss.v042.i08}
#'
#' @examples
#' data("lalonde")
#'
#' # Default: 1:1 NN PS matching w/o replacement
#'
#' m.out1 <- matchit(treat ~ age + educ + race + nodegree +
#'                    married + re74 + re75, data = lalonde)
#' m.out1
#' summary(m.out1)
#'
#' # 1:1 NN Mahalanobis distance matching w/ replacement and
#' # exact matching on married and race
#'
#' m.out2 <- matchit(treat ~ age + educ + race + nodegree +
#'                    married + re74 + re75, data = lalonde,
#'                    distance = "mahalanobis", replace = TRUE,
#'                    exact = ~ married + race)
#' m.out2
#' summary(m.out2, un = TRUE)
#'
#' # 2:1 NN Mahalanobis distance matching within caliper defined
#' # by a probit pregression PS
#'
#' m.out3 <- matchit(treat ~ age + educ + race + nodegree +
#'                    married + re74 + re75, data = lalonde,
#'                    distance = "glm", link = "probit",
#'                    mahvars = ~ age + educ + re74 + re75,
#'                    caliper = .1, ratio = 2)
#' m.out3
#' summary(m.out3, un = TRUE)
#'
#' # Optimal full PS matching for the ATE within calipers on
#' # PS, age, and educ
#' @examplesIf requireNamespace("optmatch", quietly = TRUE)
#' m.out4 <- matchit(treat ~ age + educ + race + nodegree +
#'                    married + re74 + re75, data = lalonde,
#'                    method = "full", estimand = "ATE",
#'                    caliper = c(.1, age = 2, educ = 1),
#'                    std.caliper = c(TRUE, FALSE, FALSE))
#' m.out4
#' summary(m.out4, un = TRUE)
#' @examples
#' # Subclassification on a logistic PS with 10 subclasses after
#' # discarding controls outside common support of PS
#'
#' s.out1 <- matchit(treat ~ age + educ + race + nodegree +
#'                    married + re74 + re75, data = lalonde,
#'                    method = "subclass", distance = "glm",
#'                    discard = "control", subclass = 10)
#' s.out1
#' summary(s.out1, un = TRUE)
#'
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
                    ...) {

  #Checking input format
  #data input
  mcall <- match.call()

  ## Process method
  .chk_null_or(method, chk::chk_string)
  if (length(method) == 1 && is.character(method)) {
    method <- tolower(method)
    method <- match_arg(method, c("exact", "cem", "nearest", "optimal", "full", "genetic", "subclass", "cardinality",
                                  "quick"))
    fn2 <- paste0("matchit2", method)
  }
  else if (is.null(method)) {
    fn2 <- "matchit2null"
  }
  else {
    .err("`method` must be the name of a supported matching method. See `?matchit` for allowable options")
  }

  #Process formula and data inputs
  .chk_formula(formula, sides = 2)

  tt <- terms(formula, data = data)
  treat.form <- update(tt, . ~ 0)
  treat.mf <- model.frame(treat.form, data = data, na.action = "na.pass")
  treat <- model.response(treat.mf)

  #Check and binarize treat
  treat <- check_treat(treat)
  if (length(treat) == 0) .err("the treatment cannot be `NULL`")

  names(treat) <- rownames(treat.mf)

  n.obs <- length(treat)

  #Process inputs
  ignored.inputs <- check.inputs(mcall = mcall, method = method, distance = distance, exact = exact,
                                 mahvars = mahvars, antiexact = antiexact, caliper = caliper, discard = discard,
                                 reestimate = reestimate, s.weights = s.weights, replace = replace,
                                 ratio = ratio, m.order = m.order, estimand = estimand)

  if (length(ignored.inputs) > 0) {
    for (i in ignored.inputs) assign(i, NULL)
  }

  #Process replace
  replace <- process.replace(replace, method, ...)

  #Process ratio
  ratio <- process.ratio(ratio, method, ...)

  #Process s.weights
  if (!is.null(s.weights)) {
    if (is.character(s.weights)) {
      if (is.null(data) || !is.data.frame(data)) {
        .err("if `s.weights` is specified a string, a data frame containing the named variable must be supplied to `data`")
      }
      if (!all(s.weights %in% names(data))) {
        .err("the name supplied to `s.weights` must be a variable in `data`")
      }
      s.weights.form <- reformulate(s.weights)
      s.weights <- model.frame(s.weights.form, data, na.action = "na.pass")
      if (ncol(s.weights) != 1) .err("`s.weights` can only contain one named variable")
      s.weights <- s.weights[[1]]
    }
    else if (inherits(s.weights, "formula")) {
      s.weights.form <- update(s.weights, NULL ~ .)
      s.weights <- model.frame(s.weights.form, data, na.action = "na.pass")
      if (ncol(s.weights) != 1) .err("`s.weights` can only contain one named variable")
      s.weights <- s.weights[[1]]
    }
    else if (!is.numeric(s.weights)) {
      .err("`s.weights` must be supplied as a numeric vector, string, or one-sided formula")
    }

    chk::chk_not_any_na(s.weights)
    if (length(s.weights) != n.obs) .err("`s.weights` must be the same length as the treatment vector")

    names(s.weights) <- names(treat)

  }

  #Process distance function
  is.full.mahalanobis <- FALSE
  fn1 <- NULL
  if (is.null(method) || !method %in% c("exact", "cem", "cardinality")) {
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
  if (!is.null(fn1) && fn1 == "distance2gam") {
    rlang::check_installed("mgcv")
    env <- environment(formula)
    covs.formula <- mgcv::interpret.gam(formula)$fake.formula
    environment(covs.formula) <- env
    covs.formula <- delete.response(terms(covs.formula, data = data))
  }
  else {
    covs.formula <- delete.response(terms(formula, data = data))
  }
  covs <- model.frame(covs.formula, data = data, na.action = "na.pass")
  k <- ncol(covs)
  for (i in seq_len(k)) {
    if (anyNA(covs[[i]]) || (is.numeric(covs[[i]]) && any(!is.finite(covs[[i]])))) {
      covariates.with.missingness <- names(covs)[i:k][vapply(i:k, function(j) anyNA(covs[[j]]) ||
                                                               (is.numeric(covs[[j]]) && any(!is.finite(covs[[j]]))),
                                                             logical(1L))]
      .err(paste0("Missing and non-finite values are not allowed in the covariates. Covariates with missingness or non-finite values:\n\t",
                  paste(covariates.with.missingness, collapse = ", ")), tidy = FALSE)
    }
    if (is.character(covs[[i]])) covs[[i]] <- factor(covs[[i]])
  }

  #Process exact, mahvars, and antiexact
  exactcovs <- process.variable.input(exact, data)
  exact <- attr(exactcovs, "terms")

  mahcovs <- process.variable.input(mahvars, data)
  mahvars <- attr(mahcovs, "terms")

  antiexactcovs <- process.variable.input(antiexact, data)
  antiexact <- attr(antiexactcovs, "terms")

  #Estimate distance, discard from common support, optionally re-estimate distance
  if (is.null(fn1) || is.full.mahalanobis) {
    #No distance measure
    dist.model <- distance <- link <- NULL
  }
  else if (fn1 == "distance2user") {
    dist.model <- link <- NULL
  }
  else {
    if (verbose) {
      cat("Estimating propensity scores... \n")
    }

    if (!is.null(s.weights)) {
      attr(s.weights, "in_ps") <- !distance %in% c("bart")
    }

    #Estimate distance
    if (is.null(distance.options$formula)) distance.options$formula <- formula
    if (is.null(distance.options$data)) distance.options$data <- data
    if (is.null(distance.options$verbose)) distance.options$verbose <- verbose
    if (is.null(distance.options$estimand)) distance.options$estimand <- estimand
    if (is.null(distance.options$weights) && !fn1 %in% c("distance2bart")) {
      distance.options$weights <- s.weights
    }

    if (!is.null(attr(distance, "link"))) distance.options$link <- attr(distance, "link")
    else distance.options$link <- link

    dist.out <- do.call(fn1, distance.options, quote = TRUE)

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
  if (is.null(fn1) || is.full.mahalanobis || fn1 == "distance2user") {
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
        else if (length(dim(distance.options[[i]])) == 2 && nrow(distance.options[[i]]) == n.obs) {
          distance.options[[i]] <- distance.options[[i]][!discarded,,drop = FALSE]
        }
      }
      dist.out <- do.call(fn1, distance.options, quote = TRUE)
      dist.model <- dist.out$model
      distance[!discarded] <- dist.out$distance
    }
  }

  #Process caliper
  calcovs <- NULL

  if (!is.null(caliper)) {
    caliper <- process.caliper(caliper, method, data, covs, mahcovs, distance, discarded, std.caliper)

    if (!is.null(attr(caliper, "cal.formula"))) {
      calcovs <- model.frame(attr(caliper, "cal.formula"), data, na.action = "na.pass")
      if (anyNA(calcovs)) .err("missing values are not allowed in the covariates named in `caliper`")
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

  info <- create_info(method, fn1, link, discard, replace, ratio,
                      mahalanobis = is.full.mahalanobis || !is.null(mahvars),
                      transform = attr(is.full.mahalanobis, "transform"),
                      subclass = match.out$subclass,
                      antiexact = colnames(antiexactcovs),
                      distance_is_matrix = !is.null(distance) && is.matrix(distance))

  #Create X.list for X output, removing duplicate variables
  X.list <- list(covs, exactcovs, mahcovs, calcovs, antiexactcovs)
  all.covs <- lapply(X.list, names)
  for (i in seq_along(X.list)[-1]) if (!is.null(X.list[[i]])) X.list[[i]][names(X.list[[i]]) %in% unlist(all.covs[1:(i-1)])] <- NULL
  X.list[vapply(X.list, is.null, logical(1L))] <- NULL

  ## putting all the results together
  out <- list(
    match.matrix = match.out[["match.matrix"]],
    subclass = match.out[["subclass"]],
    weights = match.out[["weights"]],
    X = do.call("cbind", X.list),
    call = mcall,
    info = info,
    estimand = estimand,
    formula = formula,
    treat = treat,
    distance = if (!is.null(distance) && !is.matrix(distance)) setNames(distance, names(treat)),
    discarded = discarded,
    s.weights = s.weights,
    exact = exact,
    mahvars = mahvars,
    caliper = caliper,
    q.cut = match.out[["q.cut"]],
    model = dist.model,
    obj = if (include.obj) match.out[["obj"]]
  )

  out[vapply(out, is.null, logical(1L))] <- NULL

  class(out) <- class(match.out)
  out
}

#' @export
#' @rdname matchit
print.matchit <- function(x, ...) {
  info <- x[["info"]]
  cal <- !is.null(x[["caliper"]])
  dis <- c("both", "control", "treat")[pmatch(info$discard, c("both", "control", "treat"), 0L)]
  disl <- length(dis) > 0
  nm <- is.null(x[["method"]])
  cat("A matchit object")
  cat(paste0("\n - method: ", info.to.method(info)))

  if (!is.null(info$distance) || info$mahalanobis) {
    cat("\n - distance: ")
    if (info$mahalanobis) {
      if (is.null(info$transform)) #mahvars used
        cat("Mahalanobis")
      else {
        cat(capwords(gsub("_", " ", info$transform, fixed = TRUE)))
      }
    }
    if (!is.null(info$distance) && !info$distance %in% matchit_distances()) {
      if (info$mahalanobis) cat(" [matching]\n             ")

      if (info$distance_is_matrix) cat("User-defined (matrix)")
      else if (info$distance != "user") cat("Propensity score")
      else if (!is.null(attr(info$distance, "custom"))) cat(attr(info$distance, "custom"))
      else cat("User-defined")

      if (cal || disl) {
        cal.ps <- "" %in% names(x[["caliper"]])
        cat(" [")
        cat(paste(c("matching", "subclassification", "caliper", "common support")[c(!nm && !info$mahalanobis && info$method != "subclass", !nm && info$method == "subclass", cal.ps, disl)], collapse = ", "))
        cat("]")
      }
      if (info$distance != "user") {
        cat("\n             - estimated with ")
        cat(info.to.distance(info))
        if (!is.null(x[["s.weights"]])) {
          if (isTRUE(attr(x[["s.weights"]], "in_ps")))
            cat("\n             - sampling weights included in estimation")
          else cat("\n             - sampling weights not included in estimation")
        }
      }
    }
  }
  if (cal) {
    cat(paste0("\n - caliper: ", paste(vapply(seq_along(x[["caliper"]]), function(z) paste0(if (names(x[["caliper"]])[z] == "") "<distance>" else names(x[["caliper"]])[z],
                                                                                            " (", format(round(x[["caliper"]][z], 3)), ")"), character(1L)),
                                       collapse = ", ")))
  }
  if (disl) {
    cat("\n - common support: ")
    if (dis == "both") cat("units from both groups")
    else if (dis == "treat") cat("treated units")
    else if (dis == "control") cat("control units")
    cat(" dropped")
  }
  cat(paste0("\n - number of obs.: ", length(x[["treat"]]), " (original)", if (!all(x[["weights"]] == 1)) paste0(", ", sum(x[["weights"]] != 0), " (matched)")))
  if (!is.null(x[["s.weights"]])) cat("\n - sampling weights: present")
  if (!is.null(x[["estimand"]])) cat(paste0("\n - target estimand: ", x[["estimand"]]))
  if (!is.null(x[["X"]])) cat(paste0("\n - covariates: ", ifelse(length(names(x[["X"]])) > 40, "too many to name", paste(names(x[["X"]]), collapse = ", "))))
  cat("\n")
  invisible(x)
}

matchit2null <- function(discarded, ...) {

  res <- list(weights = as.numeric(!discarded))
  class(res) <- "matchit"

  res
}
