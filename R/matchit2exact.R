#' Exact Matching
#' @name method_exact
#' @aliases method_exact
#' @usage NULL
#'
#' @description
#' In [matchit()], setting `method = "exact"` performs exact matching.
#' With exact matching, a complete cross of the covariates is used to form
#' subclasses defined by each combination of the covariate levels. Any subclass
#' that doesn't contain both treated and control units is discarded, leaving
#' only subclasses containing treatment and control units that are exactly
#' equal on the included covariates. The benefits of exact matching are that
#' confounding due to the covariates included is completely eliminated,
#' regardless of the functional form of the treatment or outcome models. The
#' problem is that typically many units will be discarded, sometimes
#' dramatically reducing precision and changing the target population of
#' inference. To use exact matching in combination with another matching method
#' (i.e., to exact match on some covariates and some other form of matching on
#' others), use the `exact` argument with that method.
#'
#' This page details the allowable arguments with `method = "exact"`. See
#' [matchit()] for an explanation of what each argument means in a general
#' context and how it can be specified.
#'
#' Below is how `matchit()` is used for exact matching:
#' \preformatted{
#' matchit(formula,
#'         data = NULL,
#'         method = "exact",
#'         estimand = "ATT",
#'         s.weights = NULL,
#'         verbose = FALSE,
#'         ...)
#'}
#'
#' @param formula a two-sided [formula] object containing the treatment and
#' covariates to be used in creating the subclasses defined by a full cross of
#' the covariate levels.
#' @param data a data frame containing the variables named in `formula`.
#' If not found in `data`, the variables will be sought in the
#' environment.
#' @param method set here to `"exact"`.
#' @param estimand a string containing the desired estimand. Allowable options
#' include `"ATT"`, `"ATC"`, and `"ATE"`. The estimand controls
#' how the weights are computed; see the Computing Weights section at
#' [matchit()] for details.
#' @param s.weights the variable containing sampling weights to be incorporated
#' into balance statistics. These weights do not affect the matching process.
#' @param verbose `logical`; whether information about the matching
#' process should be printed to the console.
#' @param \dots ignored.
#'
#' The arguments `distance` (and related arguments), `exact`, `mahvars`, `discard` (and related arguments), `replace`, `m.order`, `caliper` (and related arguments), and `ratio` are ignored with a warning.
#'
#' @section Outputs:
#'
#' All outputs described in [matchit()] are returned with
#' `method = "exact"` except for `match.matrix`. This is because
#' matching strata are not indexed by treated units as they are in some other
#' forms of matching. `include.obj` is ignored.
#'
#' @seealso [matchit()] for a detailed explanation of the inputs and outputs of
#' a call to `matchit()`. The `exact` argument can be used with other
#' methods to perform exact matching in combination with other matching
#' methods.
#'
#' [method_cem] for coarsened exact matching, which performs exact
#' matching on coarsened versions of the covariates.
#'
#' @references
#' In a manuscript, you don't need to cite another package when
#' using `method = "exact"` because the matching is performed completely
#' within *MatchIt*. For example, a sentence might read:
#'
#' *Exact matching was performed using the MatchIt package (Ho, Imai,
#' King, & Stuart, 2011) in R.*
#'
#' @examples
#'
#' data("lalonde")
#'
#' # Exact matching on age, race, married, and educ
#' m.out1 <- matchit(treat ~ age + race + married + educ, data = lalonde,
#'                   method = "exact")
#' m.out1
#' summary(m.out1)
#'
NULL

matchit2exact <- function(treat, covs, data, estimand = "ATT", verbose = FALSE, ...){

  if(verbose)
    cat("Exact matching... \n")

  if (length(covs) == 0) .err("covariates must be specified in the input formula to use exact matching")

  estimand <- toupper(estimand)
  estimand <- match_arg(estimand, c("ATT", "ATC", "ATE"))

  xx <- exactify(covs, names(treat))
  cc <- do.call("intersect", lapply(unique(treat), function(t) xx[treat == t]))

  if (length(cc) == 0) {
    .err("No exact matches were found")
  }

  psclass <- setNames(factor(match(xx, cc), nmax = length(cc)), names(treat))

  if (verbose) cat("Calculating matching weights... ")

  res <- list(subclass = psclass,
              weights = get_weights_from_subclass(psclass, treat, estimand))

  if (verbose) cat("Done.\n")

  class(res) <- "matchit"
  res
}
