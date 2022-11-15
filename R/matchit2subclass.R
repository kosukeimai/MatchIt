#' Subclassification
#' @rdname method_subclass
#' @aliases method_subclass
#' @usage NULL
#'
#' @description
#' In [matchit()], setting `method = "subclass"` performs
#' subclassification on the distance measure (i.e., propensity score).
#' Treatment and control units are placed into subclasses based on quantiles of
#' the propensity score in the treated group, in the control group, or overall,
#' depending on the desired estimand. Weights are computed based on the
#' proportion of treated units in each subclass. Subclassification implemented
#' here does not rely on any other package.
#'
#' This page details the allowable arguments with `method = "subclass"`.
#' See [matchit()] for an explanation of what each argument means in a general
#' context and how it can be specified.
#'
#' Below is how `matchit()` is used for subclassification:
#' \preformatted{
#' matchit(formula,
#'         data = NULL,
#'         method = "subclass",
#'         distance = "glm",
#'         link = "logit",
#'         distance.options = list(),
#'         estimand = "ATT",
#'         discard = "none",
#'         reestimate = FALSE,
#'         s.weights = NULL,
#'         verbose = FALSE,
#'         ...) }
#'
#' @param formula a two-sided [formula] object containing the treatment and
#' covariates to be used in creating the distance measure used in the
#' subclassification.
#' @param data a data frame containing the variables named in `formula`.
#' If not found in `data`, the variables will be sought in the
#' environment.
#' @param method set here to `"subclass"`.
#' @param distance the distance measure to be used. See [`distance`]
#' for allowable options. Must be a vector of distance scores or the name of a method of estimating propensity scores.
#' @param link when `distance` is specified as a string, an additional
#' argument controlling the link function used in estimating the distance
#' measure. See [`distance`] for allowable options with each option.
#' @param distance.options a named list containing additional arguments
#' supplied to the function that estimates the distance measure as determined
#' by the argument to `distance`.
#' @param estimand the target `estimand`. If `"ATT"`, the default,
#' subclasses are formed based on quantiles of the distance measure in the
#' treated group; if `"ATC"`, subclasses are formed based on quantiles of
#' the distance measure in the control group; if `"ATE"`, subclasses are
#' formed based on quantiles of the distance measure in the full sample. The
#' estimand also controls how the subclassification weights are computed; see
#' the Computing Weights section at [matchit()] for details.
#' @param discard a string containing a method for discarding units outside a
#' region of common support.
#' @param reestimate if `discard` is not `"none"`, whether to
#' re-estimate the propensity score in the remaining sample prior to
#' subclassification.
#' @param s.weights the variable containing sampling weights to be incorporated
#' into propensity score models and balance statistics.
#' @param verbose `logical`; whether information about the matching
#' process should be printed to the console.
#' @param \dots additional arguments that control the subclassification:
#' \describe{
#' \item{`subclass`}{either the number of subclasses desired
#' or a vector of quantiles used to divide the distance measure into
#' subclasses. Default is 6.}
#' \item{`min.n`}{ the minimum number of
#' units of each treatment group that are to be assigned each subclass. If the
#' distance measure is divided in such a way that fewer than `min.n` units
#' of a treatment group are assigned a given subclass, units from other
#' subclasses will be reassigned to fill the deficient subclass. Default is 1.
#' }
#' }
#'
#' The arguments `exact`, `mahvars`, `replace`, `m.order`, `caliper` (and related arguments), and `ratio` are ignored with a warning.
#'
#' @section Outputs:
#'
#' All outputs described in [matchit()] are returned with
#' `method = "subclass"` except that `match.matrix` is excluded and
#' one additional component, `q.cut`, is included, containing a vector of
#' the distance measure cutpoints used to define the subclasses. Note that when
#' `min.n > 0`, the subclass assignments may not strictly obey the
#' quantiles listed in `q.cut`. `include.obj` is ignored.
#'
#' @details
#' After subclassification, effect estimates can be computed separately in the
#' subclasses and combined, or a single marginal effect can be estimated by
#' using the weights in the full sample. When using the weights, the method is
#' sometimes referred to as marginal mean weighting through stratification
#' (MMWS; Hong, 2010) or fine stratification weighting (Desai et al., 2017).
#' The weights can be interpreted just like inverse probability weights. See `vignette("estimating-effects")` for details.
#'
#' Changing `min.n` can change the quality of the weights. Generally, a
#' low `min.w` will yield better balance because subclasses only contain
#' units with relatively similar distance values, but may yield higher variance
#' because extreme weights can occur due to there being few members of a
#' treatment group in some subclasses. When `min.n = 0`, some subclasses may fail to
#' contain units from both treatment groups, in which case all units in such subclasses
#' will be dropped.
#'
#' Note that subclassification weights can also be estimated using
#' *WeightIt*, which provides some additional methods for estimating
#' propensity scores. Where propensity score-estimation methods overlap, both
#' packages will yield the same weights.
#'
#' @seealso [matchit()] for a detailed explanation of the inputs and outputs of
#' a call to `matchit()`.
#'
#' [`method_full`] for optimal full matching, which is similar to
#' subclassification except that the number of subclasses and subclass
#' membership are chosen to optimize the within-subclass distance.
#'
#' @references In a manuscript, you don't need to cite another package when
#' using `method = "subclass"` because the subclassification is performed
#' completely within *MatchIt*. For example, a sentence might read:
#'
#' *Propensity score subclassification was performed using the MatchIt
#' package (Ho, Imai, King, & Stuart, 2011) in R.*
#'
#' It may be a good idea to cite Hong (2010) or Desai et al. (2017) if the
#' treatment effect is estimated using the subclassification weights.
#'
#' Desai, R. J., Rothman, K. J., Bateman, B. . T., Hernandez-Diaz, S., &
#' Huybrechts, K. F. (2017). A Propensity-score-based Fine Stratification
#' Approach for Confounding Adjustment When Exposure Is Infrequent:
#' Epidemiology, 28(2), 249–257. \doi{10.1097/EDE.0000000000000595}
#'
#' Hong, G. (2010). Marginal mean weighting through stratification: Adjustment
#' for selection bias in multilevel data. Journal of Educational and Behavioral
#' Statistics, 35(5), 499–531. \doi{10.3102/1076998609359785}
#'
#' @examples
#'
#' data("lalonde")
#'
#' # PS subclassification for the ATT with 7 subclasses
#' s.out1 <- matchit(treat ~ age + educ + race + nodegree +
#'                     married + re74 + re75, data = lalonde,
#'                   method = "subclass", subclass = 7)
#' s.out1
#' summary(s.out1, subclass = TRUE)
#'
#' # PS subclassification for the ATE with 10 subclasses
#' # and at least 2 units in each group per subclass
#' s.out2 <- matchit(treat ~ age + educ + race + nodegree +
#'                     married + re74 + re75, data = lalonde,
#'                   method = "subclass", subclass = 10,
#'                   estimand = "ATE", min.n = 2)
#' s.out2
#' summary(s.out2)
#'

matchit2subclass <- function(treat, distance, discarded,
                             replace = FALSE, exact = NULL,
                             estimand = "ATT", verbose = FALSE,
                             ...) {

  if(verbose)
    cat("Subclassifying... \n")

  A <- list(...)
  subclass <- A[["subclass"]]
  sub.by <- A[["sub.by"]]
  min.n <- A[["min.n"]]

  #Checks
  if (is.null(subclass)) subclass <- 6
  else if (!is.numeric(subclass) || !is.null(dim(subclass))) {
    stop("subclass must be a numeric value.", call. = FALSE)
  }
  else if (length(subclass) == 1) {
    if (round(subclass) <= 1) {
      stop("subclass must be greater than 1.",call.=FALSE)
    }
  }
  else if (!all(subclass <= 1 & subclass >= 0)) {
    stop("When specifying subclass as a vector of quantiles, all values must be between 0 and 1.",
         call. = FALSE)
  }

  if (!is.null(sub.by)) {
    sub.by.choices <- c("treat", "control", "all")
    if (!is.character(sub.by) || length(sub.by) != 1 || anyNA(pmatch(sub.by, sub.by.choices))) {
      stop("'sub.by' is deprecated and can't be converted into a proper input. Please supply an argument to 'estimand' instead.", call. = FALSE)
    }
    else {
      sub.by <- sub.by.choices[pmatch(sub.by, sub.by.choices)]
      estimand <- switch(sub.by, "treat" = "ATT", "control" = "ATC", "ATE")
      warning(sprintf("'sub.by' is deprecated and has been replaced with 'estimand'. Setting 'estimand' to \"%s\".",
                      estimand), call. = FALSE, immediate. = TRUE)
    }
  }
  else {
    estimand <- toupper(estimand)
    estimand <- match_arg(estimand, c("ATT", "ATC", "ATE"))
  }

  if (is.null(min.n)) min.n <- 1
  else if (!is.numeric(min.n) || length(min.n) != 1) {
    stop("'min.n' must be a single number.", call. = FALSE)
  }

  n.obs <- length(treat)

  ## Setting Cut Points
  if (length(subclass) == 1) {
    sprobs <- seq(0, 1, length.out = round(subclass) + 1)
  }
  else {
    sprobs <- sort(subclass)
    if (sprobs[1] != 0) sprobs <- c(0, sprobs)
    if (sprobs[length(sprobs)] != 1) sprobs <- c(sprobs, 1)
    subclass <- length(sprobs) - 1
  }

  q <- switch(estimand,
              "ATT" = quantile(distance[treat==1], probs = sprobs, na.rm = TRUE),
              "ATC" = quantile(distance[treat==0], probs = sprobs, na.rm = TRUE),
              "ATE" = quantile(distance, probs = sprobs, na.rm = TRUE))

  ## Calculating Subclasses
  psclass <- setNames(rep(NA_integer_, n.obs), names(treat))
  psclass[!discarded] <- as.integer(findInterval(distance[!discarded], q, all.inside = TRUE))

  if (length(unique(na.omit(psclass))) != subclass){
    warning("Due to discreteness in the distance measure, fewer subclasses were generated than were requested.", call.=FALSE)
  }

  if (min.n == 0) {
    ## If any subclass are missing treated or control units, set all to NA
    is.na(psclass)[!discarded & !psclass %in% intersect(psclass[!discarded & treat == 1],
                                                        psclass[!discarded & treat == 0])] <- TRUE
  }
  else if (any(table(treat, psclass) < min.n)) {
    ## If any subclasses don't have members of a treatment group, fill them
    ## by "scooting" units from nearby subclasses until each subclass has a unit
    ## from each treatment group
    psclass[!discarded] <- subclass_scoot(psclass[!discarded], treat[!discarded], distance[!discarded], min.n)
  }

  psclass <- setNames(factor(psclass, nmax = length(q)), names(treat))
  levels(psclass) <- as.character(seq_len(nlevels(psclass)))

  if (verbose) cat("Calculating matching weights... ")

  res <- list(subclass = psclass, q.cut = q,
              weights = weights.subclass(psclass, treat, estimand))

  if (verbose) cat("Done.\n")

  class(res) <- c("matchit.subclass", "matchit")
  return(res)
}
