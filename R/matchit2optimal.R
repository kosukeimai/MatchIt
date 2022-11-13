#' Optimal Pair Matching
#' @name method_optimal
#' @aliases method_optimal
#' @usage NULL
#'
#' @description
#' In [matchit()], setting `method = "optimal"` performs optimal pair
#' matching. The matching is optimal in the sense that that sum of the absolute
#' pairwise distances in the matched sample is as small as possible. The method
#' functionally relies on \pkgfun{optmatch}{fullmatch}.
#'
#' Advantages of optimal pair matching include that the matching order is not
#' required to be specified and it is less likely that extreme within-pair
#' distances will be large, unlike with nearest neighbor matching. Generally,
#' however, as a subset selection method, optimal pair matching tends to
#' perform similarly to nearest neighbor matching in that similar subsets of
#' units will be selected to be matched.
#'
#' This page details the allowable arguments with `method = "optmatch"`.
#' See [matchit()] for an explanation of what each argument means in a general
#' context and how it can be specified.
#'
#' Below is how `matchit()` is used for optimal pair matching:
#' \preformatted{
#' matchit(formula,
#'         data = NULL,
#'         method = "optimal",
#'         distance = "glm",
#'         link = "logit",
#'         distance.options = list(),
#'         estimand = "ATT",
#'         exact = NULL,
#'         mahvars = NULL,
#'         antiexact = NULL,
#'         discard = "none",
#'         reestimate = FALSE,
#'         s.weights = NULL,
#'         ratio = 1,
#'         min.controls = NULL,
#'         max.controls = NULL,
#'         verbose = FALSE,
#'         ...) }
#'
#' @param formula a two-sided [formula] object containing the treatment and
#' covariates to be used in creating the distance measure used in the matching.
#' This formula will be supplied to the functions that estimate the distance
#' measure.
#' @param data a data frame containing the variables named in `formula`.
#' If not found in `data`, the variables will be sought in the
#' environment.
#' @param method set here to `"optimal"`.
#' @param distance the distance measure to be used. See [`distance`]
#' for allowable options. Can be supplied as a distance matrix.
#' @param link when `distance` is specified as a method of estimating
#' propensity scores, an additional argument controlling the link function used
#' in estimating the distance measure. See [`distance`] for allowable
#' options with each option.
#' @param distance.options a named list containing additional arguments
#' supplied to the function that estimates the distance measure as determined
#' by the argument to `distance`.
#' @param estimand a string containing the desired estimand. Allowable options
#' include `"ATT"` and `"ATC"`. See Details.
#' @param exact for which variables exact matching should take place.
#' @param mahvars for which variables Mahalanobis distance matching should take
#' place when `distance` corresponds to a propensity score (e.g., for
#' caliper matching or to discard units for common support). If specified, the
#' distance measure will not be used in matching.
#' @param antiexact for which variables ant-exact matching should take place.
#' Anti-exact matching is processed using \pkgfun{optmatch}{antiExactMatch}.
#' @param discard a string containing a method for discarding units outside a
#' region of common support. Only allowed when `distance` is not
#' `"mahalanobis"` and not a matrix.
#' @param reestimate if `discard` is not `"none"`, whether to
#' re-estimate the propensity score in the remaining sample prior to matching.
#' @param s.weights the variable containing sampling weights to be incorporated
#' into propensity score models and balance statistics.
#' @param ratio how many control units should be matched to each treated unit
#' for k:1 matching. For variable ratio matching, see section "Variable Ratio
#' Matching" in Details below.
#' @param min.controls,max.controls for variable ratio matching, the minimum
#' and maximum number of controls units to be matched to each treated unit. See
#' section "Variable Ratio Matching" in Details below.
#' @param verbose `logical`; whether information about the matching
#' process should be printed to the console. What is printed depends on the
#' matching method. Default is `FALSE` for no printing other than
#' warnings.
#' @param \dots additional arguments passed to \pkgfun{optmatch}{fullmatch}.
#' Allowed arguments include `tol` and `solver`. See the
#' \pkgfun{optmatch}{fullmatch} documentation for details. In general, `tol`
#' should be set to a low number (e.g., `1e-7`) to get a more precise
#' solution.
#'
#' The arguments `replace`, `caliper`, and `m.order` are ignored with a warning.
#'
#' @section Outputs:
#'
#' All outputs described in [matchit()] are returned with
#' `method = "optimal"`. When `include.obj = TRUE` in the call to
#' `matchit()`, the output of the call to `optmatch::fullmatch()` will be
#' included in the output. When `exact` is specified, this will be a list
#' of such objects, one for each stratum of the `exact` variables.
#'
#' @details
#'
#' ## Mahalanobis Distance Matching
#'
#' Mahalanobis distance matching can be done one of two ways:
#'
#' 1) If no propensity score needs to be estimated, `distance` should be
#' set to `"mahalanobis"`, and Mahalanobis distance matching will occur
#' using all the variables in `formula`. Arguments to `discard` and
#' `mahvars` will be ignored. For example, to perform simple Mahalanobis
#' distance matching, the following could be run:
#'
#' \preformatted{
#' matchit(treat ~ X1 + X2, method = "nearest",
#'         distance = "mahalanobis") }
#'
#' With this code, the Mahalanobis distance is computed using `X1` and
#' `X2`, and matching occurs on this distance. The `distance`
#' component of the `matchit()` output will be empty.
#'
#' 2) If a propensity score needs to be estimated for common support with
#' `discard`, `distance` should be whatever method is used to
#' estimate the propensity score or a vector of distance measures, i.e., it
#' should not be `"mahalanobis"`. Use `mahvars` to specify the
#' variables used to create the Mahalanobis distance. For example, to perform
#' Mahalanobis after discarding units outside the common support of the
#' propensity score in both groups, the following could be run:
#'
#' \preformatted{
#' matchit(treat ~ X1 + X2 + X3, method = "nearest",
#'         distance = "glm", discard = "both",
#'         mahvars = ~ X1 + X2) }
#'
#' With this code, `X1`, `X2`, and `X3` are used to estimate the
#' propensity score (using the `"glm"` method, which by default is
#' logistic regression), which is used to identify the common support. The
#' actual matching occurs on the Mahalanobis distance computed only using
#' `X1` and `X2`, which are supplied to `mahvars`. The estimated
#' propensity scores will be included in the `distance` component of the
#' `matchit()` output.
#'
#' ## Estimand
#'
#' The `estimand` argument controls whether control units are selected to be matched with treated units
#' (`estimand = "ATT"`) or treated units are selected to be matched with
#' control units (`estimand = "ATC"`). The "focal" group (e.g., the
#' treated units for the ATT) is typically made to be the smaller treatment
#' group, and a warning will be thrown if it is not set that way unless
#' `replace = TRUE`. Setting `estimand = "ATC"` is equivalent to
#' swapping all treated and control labels for the treatment variable. When
#' `estimand = "ATC"`, the `match.matrix` component of the output
#' will have the names of the control units as the rownames and be filled with
#' the names of the matched treated units (opposite to when `estimand =
#' "ATT"`). Note that the argument supplied to `estimand` doesn't
#' necessarily correspond to the estimand actually targeted; it is merely a
#' switch to trigger which treatment group is considered "focal".
#'
#' ## Variable Ratio Matching
#'
#' `matchit()` can perform variable
#' ratio matching, which involves matching a different number of control units
#' to each treated unit. When `ratio > 1`, rather than requiring all
#' treated units to receive `ratio` matches, the arguments to
#' `max.controls` and `min.controls` can be specified to control the
#' maximum and minimum number of matches each treated unit can have.
#' `ratio` controls how many total control units will be matched: `n1 * ratio` control units will be matched, where `n1` is the number of
#' treated units, yielding the same total number of matched controls as fixed
#' ratio matching does.
#'
#' Variable ratio matching can be used with any `distance` specification.
#' `ratio` does not have to be an integer but must be greater than 1 and
#' less than `n0/n1`, where `n0` and `n1` are the number of
#' control and treated units, respectively. Setting `ratio = n0/n1`
#' performs a restricted form of full matching where all control units are
#' matched. If `min.controls` is not specified, it is set to 1 by default.
#' `min.controls` must be less than `ratio`, and `max.controls`
#' must be greater than `ratio`. See the Examples section of
#' [method_nearest()] for an example of their use, which is the same
#' as it is with optimal matching.
#'
#' @note
#' Optimal pair matching is a restricted form of optimal full matching
#' where the number of treated units in each subclass is equal to 1, whereas in
#' unrestricted full matching, multiple treated units can be assigned to the
#' same subclass. \pkgfun{optmatch}{pairmatch} is simply a wrapper for
#' \pkgfun{optmatch}{fullmatch}, which performs optimal full matching and is the
#' workhorse for [`method_full`]. In the same way, `matchit()`
#' uses `optmatch::fullmatch()` under the hood, imposing the restrictions that
#' make optimal full matching function like optimal pair matching (which is
#' simply to set `min.controls >= 1` and to pass `ratio` to the
#' `mean.controls` argument). This distinction is not important for
#' regular use but may be of interest to those examining the source code.
#'
#' The option `"optmatch_max_problem_size"` is automatically set to
#' `Inf` during the matching process, different from its default in
#' *optmatch*. This enables matching problems of any size to be run, but
#' may also let huge, infeasible problems get through and potentially take a
#' long time or crash R. See \pkgfun{optmatch}{setMaxProblemSize} for more details.
#'
#' @seealso [matchit()] for a detailed explanation of the inputs and outputs of
#' a call to `matchit()`.
#'
#' \pkgfun{optmatch}{fullmatch}, which is the workhorse.
#'
#' [`method_full`] for optimal full matching, of which optimal pair
#' matching is a special case, and which relies on similar machinery.
#'
#' @references In a manuscript, be sure to cite the following paper if using
#' `matchit()` with `method = "optimal"`:
#'
#' Hansen, B. B., & Klopfer, S. O. (2006). Optimal Full Matching and Related
#' Designs via Network Flows. Journal of Computational and Graphical
#' Statistics, 15(3), 609–627. \doi{10.1198/106186006X137047}
#'
#' For example, a sentence might read:
#'
#' *Optimal pair matching was performed using the MatchIt package (Ho,
#' Imai, King, & Stuart, 2011) in R, which calls functions from the optmatch
#' package (Hansen & Klopfer, 2006).*
#'
#' @examplesIf requireNamespace("optmatch", quietly = TRUE)
#' data("lalonde")
#'
#' # 1:1 optimal PS matching with exact matching on race
#' m.out1 <- matchit(treat ~ age + educ + race + nodegree +
#'                     married + re74 + re75, data = lalonde,
#'                   method = "optimal", exact = ~race)
#' m.out1
#' summary(m.out1)
#'
#' #2:1 optimal Mahalanobis distance matching
#' m.out2 <- matchit(treat ~ age + educ + race + nodegree +
#'                     married + re74 + re75, data = lalonde,
#'                   method = "optimal", distance = "mahalanobis",
#'                   ratio = 2)
#' m.out2
#' summary(m.out2, un = FALSE)
NULL

matchit2optimal <- function(treat, formula, data, distance, discarded,
                            ratio = 1, s.weights = NULL, caliper = NULL,
                            mahvars = NULL, exact = NULL,
                            estimand = "ATT", verbose = FALSE,
                            is.full.mahalanobis,  antiexact = NULL, ...) {

  check.package("optmatch")

  if (verbose) cat("Optimal matching... \n")

  A <- list(...)
  pm.args <- c("tol", "solver")
  A[!names(A) %in% pm.args] <- NULL

  #Set max problem size to Inf and return to original value after match
  omps <- getOption("optmatch_max_problem_size")
  on.exit(options(optmatch_max_problem_size = omps))
  options(optmatch_max_problem_size = Inf)

  estimand <- toupper(estimand)
  estimand <- match_arg(estimand, c("ATT", "ATC"))
  if (estimand == "ATC") {
    tc <- c("control", "treated")
    focal <- 0
  }
  else {
    tc <- c("treated", "control")
    focal <- 1
  }

  treat_ <- setNames(as.integer(treat[!discarded] == focal), names(treat)[!discarded])

  # if (!is.null(data)) data <- data[!discarded,]

  if (is.full.mahalanobis) {
    if (length(attr(terms(formula, data = data), "term.labels")) == 0) {
      stop(sprintf("Covariates must be specified in the input formula when distance = \"%s\".",
                   attr(is.full.mahalanobis, "transform")), call. = FALSE)
    }
    mahvars <- formula
  }

  if (!is.null(caliper)) {
    warning("Calipers are currently not compatible with method = \"optimal\" and will be ignored.", call. = FALSE, immediate. = TRUE)
    caliper <- NULL
  }

  min.controls <- attr(ratio, "min.controls")
  max.controls <- attr(ratio, "max.controls")

  if (is.null(max.controls)) {
    min.controls <- max.controls <- ratio
  }

  #Exact matching strata
  if (!is.null(exact)) {
    ex <- factor(exactify(model.frame(exact, data = data),
                          sep = ", ", include_vars = TRUE)[!discarded])

    cc <- intersect(ex[treat_==1], ex[treat_==0])
    if (length(cc) == 0) stop("No matches were found.", call. = FALSE)

    e_ratios <- vapply(levels(ex), function(e) sum(treat_[ex == e] == 0)/sum(treat_[ex == e] == 1), numeric(1L))

    if (any(e_ratios < 1)) {
      warning(sprintf("Fewer %s units than %s units in some 'exact' strata; not all %s units will get a match.",
                      tc[2], tc[1], tc[1]), immediate. = TRUE, call. = FALSE)
    }
    if (ratio > 1 && any(e_ratios < ratio)) {
      if (ratio == max.controls)
        warning(sprintf("Not all %s units will get %s matches.",
                        tc[1], ratio), immediate. = TRUE, call. = FALSE)
      else
        warning(sprintf("Not enough %s units for an average of %s matches per %s unit in all 'exact' strata.",
                        tc[2], ratio, tc[1]), immediate. = TRUE, call. = FALSE)
    }
  }
  else {
    ex <- factor(rep("_", length(treat_)), levels = "_")
    e_ratios <- setNames(sum(treat_ == 0)/sum(treat_ == 1), levels(ex))

    if (e_ratios < 1) {
      warning(sprintf("Fewer %s units than %s units; not all %s units will get a match.",
                      tc[2], tc[1], tc[1]), immediate. = TRUE, call. = FALSE)
    }
    else if (e_ratios < ratio) {
      if (ratio == max.controls)
        warning(sprintf("Not all %s units will get %s matches.",
                        tc[1], ratio), immediate. = TRUE, call. = FALSE)
      else
        warning(sprintf("Not enough %s units for an average of %s matches per %s unit.",
                        tc[2], ratio, tc[1]), immediate. = TRUE, call. = FALSE)
    }
  }

  #Create distance matrix; note that Mahalanobis distance computed using entire
  #sample (minus discarded), like method2nearest, as opposed to within exact strata, like optmatch.
  if (!is.null(mahvars)) {
    transform <- if (is.full.mahalanobis) attr(is.full.mahalanobis, "transform") else "mahalanobis"
    mahcovs <- transform_covariates(mahvars, data = data, method = transform,
                                    s.weights = s.weights, treat = treat,
                                    discarded = discarded)
    mo <- eucdist_internal(mahcovs, treat)
  }
  else if (is.matrix(distance)) {
    mo <- distance
  }
  else {
    mo <- eucdist_internal(setNames(distance, names(treat)), treat)
  }

  #Transpose distance mat as needed
  if (focal == 0) mo <- t(mo)

  #Remove discarded units from distance mat
  mo <- mo[!discarded[treat == focal], !discarded[treat != focal], drop = FALSE]
  dimnames(mo) <- list(names(treat_)[treat_ == 1], names(treat_)[treat_ == 0])

  mo <- optmatch::match_on(mo, data = data[!discarded,, drop = FALSE])
  mo <- optmatch::as.InfinitySparseMatrix(mo)

  #Process antiexact
  if (!is.null(antiexact)) {
    antiexactcovs <- model.frame(antiexact, data)
    for (i in seq_len(ncol(antiexactcovs))) {
      mo <- mo + optmatch::antiExactMatch(antiexactcovs[[i]][!discarded], z = treat_)
    }
  }

  #Initialize pair membership; must include names
  pair <- setNames(rep(NA_character_, length(treat)), names(treat))
  p <- setNames(vector("list", nlevels(ex)), levels(ex))

  t_df <- data.frame(treat)

  for (e in levels(ex)) {
    if (nlevels(ex) > 1) {
      mo_ <- mo[ex[treat_==1] == e, ex[treat_==0] == e]
    }
    else mo_ <- mo

    if (any(dim(mo_) == 0)) next
    else if (all(dim(mo_) == 1)) {
      pair[ex == e] <- paste(1, e, sep = "|")
      next
    }

    #Process ratio, etc., when available ratio in exact matching categories
    #(e_ratio) differs from requested ratio
    if (e_ratios[e] < 1) {
      #Switch treatment and control labels; unmatched treated units are dropped
      ratio_ <- min.controls_ <- max.controls_ <- 1
      mo_ <- t(mo_)
    }
    else if (e_ratios[e] < ratio) {
      #Lower ratio and min.controls.
      ratio_ <- e_ratios[e]
      min.controls_ <- min(min.controls, floor(e_ratios[e]))
      max.controls_ <- max.controls
    }
    else {
      ratio_ <- ratio
      min.controls_ <- min.controls
      max.controls_ <- max.controls
    }

    withCallingHandlers({
      p[[e]] <- do.call(optmatch::fullmatch,
                        c(list(mo_,
                               mean.controls = ratio_,
                               min.controls = min.controls_,
                               max.controls = max.controls_,
                               data = t_df), #just to get rownames; not actually used in matching
                          A))
    },
    warning = function(w) {
      warning(paste0("(from optmatch) ", conditionMessage(w)), call. = FALSE, immediate. = TRUE)
      invokeRestart("muffleWarning")
    },
    error = function(e1) {
      stop(paste0("(from optmatch) ", conditionMessage(e1)), call. = FALSE)
    })

    pair[names(p[[e]])[!is.na(p[[e]])]] <- paste(as.character(p[[e]][!is.na(p[[e]])]), e, sep = "|")
  }

  if (all(is.na(pair))) stop("No matches were found.", call. = FALSE)
  if (length(p) == 1) p <- p[[1]]

  psclass <- factor(pair)
  levels(psclass) <- seq_len(nlevels(psclass))
  names(psclass) <- names(treat)

  mm <- nummm2charmm(subclass2mmC(psclass, treat, focal), treat)

  if (verbose) cat("Calculating matching weights... ")

  ## calculate weights and return the results
  res <- list(match.matrix = mm,
              subclass = psclass,
              weights = weights.subclass(psclass, treat, estimand),
              obj = p)

  if (verbose) cat("Done.\n")

  class(res) <- "matchit"
  return(res)
}
