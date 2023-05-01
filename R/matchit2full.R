#' Optimal Full Matching
#' @name method_full
#' @aliases method_full
#' @usage NULL
#'
#' @description
#' In [matchit()], setting `method = "full"` performs optimal full
#' matching, which is a form of subclassification wherein all units, both
#' treatment and control (i.e., the "full" sample), are assigned to a subclass
#' and receive at least one match. The matching is optimal in the sense that
#' that sum of the absolute distances between the treated and control units in
#' each subclass is as small as possible. The method relies on and is a wrapper
#' for \pkgfun{optmatch}{fullmatch}.
#'
#' Advantages of optimal full matching include that the matching order is not
#' required to be specified, units do not need to be discarded, and it is less
#' likely that extreme within-subclass distances will be large, unlike with
#' standard subclassification. The primary output of full matching is a set of
#' matching weights that can be applied to the matched sample; in this way,
#' full matching can be seen as a robust alternative to propensity score
#' weighting, robust in the sense that the propensity score model does not need
#' to be correct to estimate the treatment effect without bias. Note: with large samples, the optimization may fail or run very slowly; one can try using [`method = "quick"`][method_quick] instead, which also performs full matching but can be much faster.
#'
#' This page details the allowable arguments with `method = "full"`.
#' See [matchit()] for an explanation of what each argument means in a general
#' context and how it can be specified.
#'
#' Below is how `matchit()` is used for optimal full matching:
#' \preformatted{
#' matchit(formula,
#'         data = NULL,
#'         method = "full",
#'         distance = "glm",
#'         link = "logit",
#'         distance.options = list(),
#'         estimand = "ATT",
#'         exact = NULL,
#'         mahvars = NULL,
#'         anitexact = NULL,
#'         discard = "none",
#'         reestimate = FALSE,
#'         s.weights = NULL,
#'         caliper = NULL,
#'         std.caliper = TRUE,
#'         verbose = FALSE,
#'         ...)
#' }
#'
#' @param formula a two-sided [formula] object containing the treatment and
#' covariates to be used in creating the distance measure used in the matching.
#' This formula will be supplied to the functions that estimate the distance
#' measure.
#' @param data a data frame containing the variables named in `formula`.
#' If not found in `data`, the variables will be sought in the
#' environment.
#' @param method set here to `"full"`.
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
#' include `"ATT"`, `"ATC"`, and `"ATE"`. The estimand controls
#' how the weights are computed; see the Computing Weights section at
#' [matchit()] for details.
#' @param exact for which variables exact matching should take place.
#' @param mahvars for which variables Mahalanobis distance matching should take
#' place when `distance` corresponds to a propensity score (e.g., for
#' caliper matching or to discard units for common support). If specified, the
#' distance measure will not be used in matching.
#' @param antiexact for which variables ant-exact matching should take place.
#' Anti-exact matching is processed using \pkgfun{optmatch}{antiExactMatch}.
#' @param discard a string containing a method for discarding units outside a
#' region of common support. Only allowed when `distance` corresponds to a
#' propensity score.
#' @param reestimate if `discard` is not `"none"`, whether to
#' re-estimate the propensity score in the remaining sample prior to matching.
#' @param s.weights the variable containing sampling weights to be incorporated
#' into propensity score models and balance statistics.
#' @param caliper the width(s) of the caliper(s) used for caliper matching.
#' Calipers are processed by \pkgfun{optmatch}{caliper}. See Notes and Examples.
#' @param std.caliper `logical`; when calipers are specified, whether they
#' are in standard deviation units (`TRUE`) or raw units (`FALSE`).
#' @param verbose `logical`; whether information about the matching
#' process should be printed to the console.
#' @param \dots additional arguments passed to \pkgfun{optmatch}{fullmatch}.
#' Allowed arguments include `min.controls`, `max.controls`,
#' `omit.fraction`, `mean.controls`, `tol`, and `solver`.
#' See the \pkgfun{optmatch}{fullmatch} documentation for details. In general,
#' `tol` should be set to a low number (e.g., `1e-7`) to get a more
#' precise solution.
#'
#' The arguments `replace`, `m.order`, and `ratio` are ignored with a warning.
#'
#' @section Outputs: All outputs described in [matchit()] are returned with
#' `method = "full"` except for `match.matrix`. This is because
#' matching strata are not indexed by treated units as they are in some other
#' forms of matching. When `include.obj = TRUE` in the call to
#' `matchit()`, the output of the call to \pkgfun{optmatch}{fullmatch} will be
#' included in the output. When `exact` is specified, this will be a list
#' of such objects, one for each stratum of the `exact` variables.
#'
#' @details
#' ## Mahalanobis Distance Matching
#'
#' Mahalanobis distance matching can be done one of two ways:
#'
#' \enumerate{
#' \item{
#' If no propensity score needs to be estimated, `distance` should be
#' set to `"mahalanobis"`, and Mahalanobis distance matching will occur
#' using all the variables in `formula`. Arguments to `discard` and
#' `mahvars` will be ignored, and a caliper can only be placed on named
#' variables. For example, to perform simple Mahalanobis distance matching, the
#' following could be run:
#'
#' \preformatted{
#' matchit(treat ~ X1 + X2, method = "nearest",
#'         distance = "mahalanobis") }
#'
#' With this code, the Mahalanobis distance is computed using `X1` and
#' `X2`, and matching occurs on this distance. The `distance`
#' component of the `matchit()` output will be empty.
#' }
#' \item{
#' If a propensity score needs to be estimated for any reason, e.g., for
#' common support with `discard` or for creating a caliper,
#' `distance` should be whatever method is used to estimate the propensity
#' score or a vector of distance measures, i.e., it should not be
#' `"mahalanobis"`. Use `mahvars` to specify the variables used to
#' create the Mahalanobis distance. For example, to perform Mahalanobis within
#' a propensity score caliper, the following could be run:
#'
#' \preformatted{
#' matchit(treat ~ X1 + X2 + X3, method = "nearest",
#'         distance =  "glm", caliper = .25,
#'         mahvars = ~ X1 + X2) }
#'
#' With this code, `X1`, `X2`, and `X3` are used to estimate the
#' propensity score (using the `"glm"` method, which by default is
#' logistic regression), which is used to create a matching caliper. The actual
#' matching occurs on the Mahalanobis distance computed only using `X1`
#' and `X2`, which are supplied to `mahvars`. Units whose propensity
#' score difference is larger than the caliper will not be paired, and some
#' treated units may therefore not receive a match. The estimated propensity
#' scores will be included in the `distance` component of the
#' `matchit()` output. See Examples.
#' }
#' }
#'
#' @note Calipers can only be used when `min.controls` is left at its
#' default.
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
#' [`method_optimal`] for optimal pair matching, which is a special
#' case of optimal full matching, and which relies on similar machinery.
#' Results from `method = "optimal"` can be replicated with `method = "full"` by setting `min.controls`, `max.controls`, and
#' `mean.controls` to the desired `ratio`.
#'
#' [`method_quick`] for fast generalized quick matching, which is very similar to optimal full matching but can be dramatically faster at the expense of optimality and is less customizable.
#'
#' @references In a manuscript, be sure to cite the following paper if using
#' `matchit()` with `method = "full"`:
#'
#' Hansen, B. B., & Klopfer, S. O. (2006). Optimal Full Matching and Related
#' Designs via Network Flows. *Journal of Computational and Graphical Statistics*,
#' 15(3), 609–627. \doi{10.1198/106186006X137047}
#'
#' For example, a sentence might read:
#'
#' *Optimal full matching was performed using the MatchIt package (Ho,
#' Imai, King, & Stuart, 2011) in R, which calls functions from the optmatch
#' package (Hansen & Klopfer, 2006).*
#'
#' Theory is also developed in the following article:
#'
#' Hansen, B. B. (2004). Full Matching in an Observational Study of Coaching
#' for the SAT. Journal of the American Statistical Association, 99(467),
#' 609–618. \doi{10.1198/016214504000000647}
#'
#' @examplesIf requireNamespace("optmatch", quietly = TRUE)
#' data("lalonde")
#'
#' # Optimal full PS matching
#' m.out1 <- matchit(treat ~ age + educ + race + nodegree +
#'                     married + re74 + re75, data = lalonde,
#'                   method = "full")
#' m.out1
#' summary(m.out1)
#'
#' # Optimal full Mahalanobis distance matching within a PS caliper
#' m.out2 <- matchit(treat ~ age + educ + race + nodegree +
#'                     married + re74 + re75, data = lalonde,
#'                   method = "full", caliper = .01,
#'                   mahvars = ~ age + educ + re74 + re75)
#' m.out2
#' summary(m.out2, un = FALSE)
#'
#' # Optimal full Mahalanobis distance matching within calipers
#' # of 500 on re74 and re75
#' m.out3 <- matchit(treat ~ age + educ + re74 + re75,
#'                   data = lalonde, distance = "mahalanobis",
#'                   method = "full",
#'                   caliper = c(re74 = 500, re75 = 500),
#'                   std.caliper = FALSE)
#' m.out3
#' summary(m.out3, addlvariables = ~race + nodegree + married,
#'         data = lalonde, un = FALSE)
NULL

matchit2full <- function(treat, formula, data, distance, discarded,
                         ratio = NULL, s.weights = NULL, #min.controls and max.controls in attrs of replace
                         caliper = NULL, mahvars = NULL, exact = NULL,
                         estimand = "ATT", verbose = FALSE,
                         is.full.mahalanobis, antiexact = NULL, ...) {

  rlang::check_installed("optmatch")

  if (verbose) cat("Full matching... \n")

  A <- list(...)

  fm.args <- c("omit.fraction", "mean.controls", "tol", "solver")
  A[!names(A) %in% fm.args] <- NULL

  #Set max problem size to Inf and return to original value after match
  omps <- getOption("optmatch_max_problem_size")
  on.exit(options(optmatch_max_problem_size = omps))
  options(optmatch_max_problem_size = Inf)

  estimand <- toupper(estimand)
  estimand <- match_arg(estimand, c("ATT", "ATC", "ATE"))
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
      .err(sprintf("covariates must be specified in the input formula when `distance = \"%s\"`",
                   attr(is.full.mahalanobis, "transform")))
    }
    mahvars <- formula
  }

  min.controls <- attr(ratio, "min.controls")
  max.controls <- attr(ratio, "max.controls")

  #Exact matching strata
  if (!is.null(exact)) {
    ex <- factor(exactify(model.frame(exact, data = data),
                          sep = ", ", include_vars = TRUE)[!discarded])
    cc <- intersect(as.integer(ex)[treat_==1], as.integer(ex)[treat_==0])
    if (length(cc) == 0) .err("No matches were found")
  }
  else {
    ex <- factor(rep("_", length(treat_)), levels = "_")
    cc <- 1
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

  #Process caliper
  if (!is.null(caliper)) {
    if (min.controls != 0) {
      .err("calipers cannot be used with `method = \"full\"` when `min.controls` is specified")
    }

    if (any(names(caliper) != "")) {
      cov.cals <- setdiff(names(caliper), "")
      calcovs <- get.covs.matrix(reformulate(cov.cals, intercept = FALSE), data = data)
    }
    for (i in seq_along(caliper)) {
      if (names(caliper)[i] != "") {
        mo_cal <- optmatch::match_on(setNames(calcovs[!discarded, names(caliper)[i]], names(treat_)), z = treat_)
      }
      else if (is.null(mahvars) || is.matrix(distance)) {
        mo_cal <- mo
      }
      else {
        mo_cal <- optmatch::match_on(setNames(distance[!discarded], names(treat_)), z = treat_)
      }

      mo <- mo + optmatch::caliper(mo_cal, caliper[i])
    }
    rm(mo_cal)
  }

  #Initialize pair membership; must include names
  pair <- setNames(rep(NA_character_, length(treat)), names(treat))
  p <- setNames(vector("list", nlevels(ex)), levels(ex))

  t_df <- data.frame(treat)

  for (e in levels(ex)[cc]) {
    if (nlevels(ex) > 1) {
      if (verbose) {
        cat(sprintf("Matching subgroup %s/%s: %s...\n",
                    match(e, levels(ex)[cc]), length(cc), e))
      }
      mo_ <- mo[ex[treat_==1] == e, ex[treat_==0] == e]
    }
    else mo_ <- mo

    if (any(dim(mo_) == 0) || !any(is.finite(mo_))) next
    else if (all(dim(mo_) == 1) && all(is.finite(mo_))) {
      pair[ex == e] <- paste(1, e, sep = "|")
      next
    }

    matchit_try({
      p[[e]] <- do.call(optmatch::fullmatch,
                        c(list(mo_,
                               min.controls = min.controls,
                               max.controls = max.controls,
                               data = t_df), #just to get rownames; not actually used in matching
                          A))
    }, from = "optmatch")

    pair[names(p[[e]])[!is.na(p[[e]])]] <- paste(as.character(p[[e]][!is.na(p[[e]])]), e, sep = "|")
  }

  if (all(is.na(pair))) .err("No matches were found")
  if (length(p) == 1) p <- p[[1]]

  psclass <- factor(pair)
  levels(psclass) <- seq_len(nlevels(psclass))
  names(psclass) <- names(treat)

  #No match.matrix because treated units don't index matched strata (i.e., more than one
  #treated unit can be in the same stratum). Stratum information is contained in subclass.

  if (verbose) cat("Calculating matching weights... ")

  res <- list(subclass = psclass,
              weights = get_weights_from_subclass(psclass, treat, estimand),
              obj = p)

  if (verbose) cat("Done.\n")

  class(res) <- c("matchit")
  res
}
