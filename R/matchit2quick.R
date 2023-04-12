#' Fast Generalized Full Matching
#' @name method_quick
#' @aliases method_quick
#' @usage NULL
#'
#' @description
#' In [matchit()], setting `method = "quick"` performs generalized full
#' matching, which is a form of subclassification wherein all units, both
#' treatment and control (i.e., the "full" sample), are assigned to a subclass
#' and receive at least one match. It uses an algorithm that is extremely fast
#' compared to optimal full matching, which is why it is labelled as "quick", at the
#' expense of true optimality. The method is described in Sävje, Higgins, & Sekhon (2021). The method relies on and is a wrapper
#' for \pkgfun{quickmatch}{quickmatch}.
#'
#' Advantages of generalized full matching include that the matching order is not
#' required to be specified, units do not need to be discarded, and it is less
#' likely that extreme within-subclass distances will be large, unlike with
#' standard subclassification. The primary output of generalized full matching is a set of
#' matching weights that can be applied to the matched sample; in this way,
#' generalized full matching can be seen as a robust alternative to propensity score
#' weighting, robust in the sense that the propensity score model does not need
#' to be correct to estimate the treatment effect without bias.
#'
#' This page details the allowable arguments with `method = "quick"`.
#' See [matchit()] for an explanation of what each argument means in a general
#' context and how it can be specified.
#'
#' Below is how `matchit()` is used for generalized full matching:
#' \preformatted{
#' matchit(formula,
#'         data = NULL,
#'         method = "quick",
#'         distance = "glm",
#'         link = "logit",
#'         distance.options = list(),
#'         estimand = "ATT",
#'         exact = NULL,
#'         mahvars = NULL,
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
#' @param method set here to `"quick"`.
#' @param distance the distance measure to be used. See [`distance`]
#' for allowable options. Cannot be supplied as a matrix.
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
#' place when `distance` corresponds to a propensity score (e.g., to discard units for common support). If specified, the
#' distance measure will not be used in matching.
#' @param discard a string containing a method for discarding units outside a
#' region of common support. Only allowed when `distance` corresponds to a
#' propensity score.
#' @param reestimate if `discard` is not `"none"`, whether to
#' re-estimate the propensity score in the remaining sample prior to matching.
#' @param s.weights the variable containing sampling weights to be incorporated
#' into propensity score models and balance statistics.
#' @param caliper the width of the caliper used for caliper matching. A caliper can only be placed on the propensity score.
#' @param std.caliper `logical`; when a caliper is specified, whether it
#' is in standard deviation units (`TRUE`) or raw units (`FALSE`).
#' @param verbose `logical`; whether information about the matching
#' process should be printed to the console.
#' @param \dots additional arguments passed to \pkgfun{quickmatch}{quickmatch}. Allowed arguments include `treatment_constraints`, `size_constraint`, `target`, and other arguments passed to `scclust::sc_clustering()` (see \pkgfun{quickmatch}{quickmatch} for details). In particular, changing `seed_method` from its default can improve performance.
#' No arguments will be passed to `distances::distances()`.
#'
#' The arguments `replace`, `ratio`, `min.controls`, `max.controls`, `m.order`, and `antiexact` are ignored with a warning.
#'
#' @section Outputs:
#'
#' All outputs described in [matchit()] are returned with
#' `method = "quick"` except for `match.matrix`. This is because
#' matching strata are not indexed by treated units as they are in some other
#' forms of matching. When `include.obj = TRUE` in the call to
#' `matchit()`, the output of the call to \pkgfun{quickmatch}{quickmatch} will be
#' included in the output. When `exact` is specified, this will be a list
#' of such objects, one for each stratum of the `exact` variables.
#'
#' @details
#'
#' Generalized full matching is similar to optimal full matching, but has some additional flexibility that can be controlled by some of the extra arguments available. By default, `method = "quick"` performs a standard full match in which all units are matched (unless restricted by the caliper) and assigned to a subclass. Each subclass could contain multiple units from each treatment group. The subclasses are chosen to minimize the largest within-subclass distance between units (including between units of the same treatment group). Notably, generalized full matching requires less memory and can run much faster than optimal full matching and optimal pair matching and, in some cases, even than nearest neighbor matching, and it can be used with huge datasets (e.g., in the millions) while running in under a minute.
#'
#'
#' @references In a manuscript, be sure to cite the *quickmatch* package if using
#' `matchit()` with `method = "quick"`:
#'
#' Sävje, F., Sekhon, J., & Higgins, M. (2018). quickmatch: Quick generalized full matching. \url{https://CRAN.R-project.org/package=quickmatch}
#'
#' For example, a sentence might read:
#'
#' *Generalized full matching was performed using the MatchIt package (Ho,
#' Imai, King, & Stuart, 2011) in R, which calls functions from the quickmatch
#' package (Savje, Sekhon, & Higgins, 2018).*
#'
#' You should also cite the following paper, which develops and describes the method:
#'
#' Sävje, F., Higgins, M. J., & Sekhon, J. S. (2021). Generalized Full Matching. *Political Analysis*, 29(4), 423–447. \doi{10.1017/pan.2020.32}
#'
#' @seealso [matchit()] for a detailed explanation of the inputs and outputs of
#' a call to `matchit()`.
#'
#' \pkgfun{quickmatch}{quickmatch}, which is the workhorse.
#'
#' [`method_full`] for optimal full matching, which is nearly the same but offers more customizability and more optimal solutions at the cost of speed.
#'
#' @examplesIf requireNamespace("quickmatch", quietly = TRUE)
#' data("lalonde")
#'
#' # Generalize full PS matching
#' m.out1 <- matchit(treat ~ age + educ + race + nodegree +
#'                     married + re74 + re75, data = lalonde,
#'                   method = "quick")
#' m.out1
#' summary(m.out1)
NULL

matchit2quick <- function(treat, formula, data, distance, discarded,
                          s.weights = NULL,
                          caliper = NULL, mahvars = NULL, exact = NULL,
                          estimand = "ATT", verbose = FALSE,
                          is.full.mahalanobis, ...) {

  rlang::check_installed("quickmatch")

  if (verbose) cat("Generalized full matching... \n")

  A <- list(...)

  distances.args <- c("data", "id_variable", "dist_variables", "normalize", "weights")
  A[names(A) %in% distances.args] <- NULL

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

  treat_ <- treat[!discarded]
  # treat_ <- setNames(as.integer(treat[!discarded] == focal), names(treat)[!discarded])

  if (is.full.mahalanobis) {
    if (length(attr(terms(formula, data = data), "term.labels")) == 0) {
      .err(sprintf("covariates must be specified in the input formula when `distance = \"%s\"`",
                   attr(is.full.mahalanobis, "transform")))
    }
    mahvars <- formula
  }

  #Exact matching strata
  if (!is.null(exact)) {
    ex <- factor(exactify(model.frame(exact, data = data),
                          sep = ", ", include_vars = TRUE)[!discarded])
    cc <- intersect(as.integer(ex)[treat_==1], as.integer(ex)[treat_==0])
    if (length(cc) == 0) .err("no matches were found")

  }
  else {
    ex <- factor(rep("_", length(treat_)), levels = "_")
    cc <- 1
  }

  #Create distance matrix; note that Mahalanobis distance computed using entire
  #sample (minus discarded), like method2nearest, as opposed to within exact strata, like optmatch.
  if (!is.null(mahvars)) {
    transform <- if (is.full.mahalanobis) attr(is.full.mahalanobis, "transform") else "mahalanobis"
    distcovs <- transform_covariates(mahvars, data = data, method = transform,
                                     s.weights = s.weights, treat = treat,
                                     discarded = discarded)
  }
  else {
    distcovs <- as.matrix(distance)
  }


  #Remove discarded units from distance mat
  distcovs <- distcovs[!discarded, , drop = FALSE]
  rownames(distcovs) <- names(treat_)

  #Process caliper
  if (!is.null(caliper)) {
    if (!is.null(mahvars)) {
      .err("with `method = \"quick\"`, a caliper can only be used when `distance` is a propensity score or vector and `mahvars` is not specified")
    }
    if (length(caliper) > 1 || !identical(names(caliper), "")) {
      .err("with `method = \"quick\"`, calipers cannot be placed on covariates")
    }
  }

  #Initialize pair membership; must include names
  pair <- setNames(rep(NA_character_, length(treat)), names(treat))
  p <- setNames(vector("list", nlevels(ex)), levels(ex))

  for (e in levels(ex)[cc]) {

    distcovs_ <- distcovs[ex == e,, drop = FALSE]

    matchit_try({
      p[[e]] <- do.call(quickmatch::quickmatch,
                        c(list(distcovs_,
                               treatments = treat_[ex == e],
                               caliper = caliper),
                          A))
    }, from = "quickmatch")

    pair[which(ex == e)[!is.na(p[[e]])]] <- paste(as.character(p[[e]][!is.na(p[[e]])]), e, sep = "|")
  }

  if (all(is.na(pair))) .err("no matches were found")
  if (length(p) == 1) p <- p[[1]]

  psclass <- factor(pair)
  levels(psclass) <- seq_len(nlevels(psclass))
  names(psclass) <- names(treat)

  #No match.matrix because treated units don't index matched strata (i.e., more than one
  #treated unit can be in the same stratum). Stratum information is contained in subclass.

  if (verbose) cat("Calculating matching weights... ")

  res <- list(subclass = psclass,
              weights = weights.subclass(psclass, treat, estimand),
              obj = p)

  if (verbose) cat("Done.\n")

  class(res) <- c("matchit")
  res
}