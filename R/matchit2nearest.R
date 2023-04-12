#' Nearest Neighbor Matching
#' @name method_nearest
#' @aliases method_nearest
#' @usage NULL
#'
#' @description
#' In [matchit()], setting `method = "nearest"` performs greedy nearest
#' neighbor matching. A distance is computed between each treated unit and each
#' control unit, and, one by one, each treated unit is assigned a control unit
#' as a match. The matching is "greedy" in the sense that there is no action
#' taken to optimize an overall criterion; each match is selected without
#' considering the other matches that may occur subsequently.
#'
#' This page details the allowable arguments with `method = "nearest"`.
#' See [matchit()] for an explanation of what each argument means in a general
#' context and how it can be specified.
#'
#' Below is how `matchit()` is used for nearest neighbor matching:
#' \preformatted{
#' matchit(formula,
#'         data = NULL,
#'         method = "nearest",
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
#'         replace = TRUE,
#'         m.order = NULL,
#'         caliper = NULL,
#'         ratio = 1,
#'         min.controls = NULL,
#'         max.controls = NULL,
#'         verbose = FALSE,
#'         ...) }
#'
#' @param formula a two-sided [formula] object containing the treatment and
#' covariates to be used in creating the distance measure used in the matching.
#' @param data a data frame containing the variables named in `formula`.
#' If not found in `data`, the variables will be sought in the
#' environment.
#' @param method set here to `"nearest"`.
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
#' @param discard a string containing a method for discarding units outside a
#' region of common support. Only allowed when `distance` corresponds to a
#' propensity score.
#' @param reestimate if `discard` is not `"none"`, whether to
#' re-estimate the propensity score in the remaining sample prior to matching.
#' @param s.weights the variable containing sampling weights to be incorporated
#' into propensity score models and balance statistics.
#' @param replace whether matching should be done with replacement.
#' @param m.order the order that the matching takes place. The default is
#' `"largest"` when `distance` corresponds to a propensity score and
#' `"data"` otherwise. See [matchit()] for allowable options.
#' @param caliper the width(s) of the caliper(s) used for caliper matching. See
#' Details and Examples.
#' @param std.caliper `logical`; when calipers are specified, whether they
#' are in standard deviation units (`TRUE`) or raw units (`FALSE`).
#' @param ratio how many control units should be matched to each treated unit
#' for k:1 matching. For variable ratio matching, see section "Variable Ratio
#' Matching" in Details below.
#' @param min.controls,max.controls for variable ratio matching, the minimum
#' and maximum number of controls units to be matched to each treated unit. See
#' section "Variable Ratio Matching" in Details below.
#' @param verbose `logical`; whether information about the matching
#' process should be printed to the console. When `TRUE`, a progress bar
#' implemented using *RcppProgress* will be displayed.
#' @param \dots additional arguments that control the matching specification:
#' \describe{
#' \item{`reuse.max`}{ `numeric`; the maximum number of
#' times each control can be used as a match. Setting `reuse.max = 1`
#' corresponds to matching without replacement (i.e., `replace = FALSE`),
#' and setting `reuse.max = Inf` corresponds to traditional matching with
#' replacement (i.e., `replace = TRUE`) with no limit on the number of
#' times each control unit can be matched. Other values restrict the number of
#' times each control can be matched when matching with replacement.
#' `replace` is ignored when `reuse.max` is specified.  }
#' \item{`unit.id`}{ one or more variables containing a unit ID for each
#' observation, i.e., in case multiple observations correspond to the same
#' unit. Once a control observation has been matched, no other observation with
#' the same unit ID can be used as matches. This ensures each control unit is
#' used only once even if it has multiple observations associated with it.
#' Omitting this argument is the same as giving each observation a unique ID.
#' Ignored when `replace = TRUE`.  }
#' }
#'
#' @note Sometimes an error will be produced by *Rcpp* along the lines of
#' `"function 'Rcpp_precious_remove' not provided by package 'Rcpp'"`. It
#' is not immediately clear why this happens, though
#' [this](https://lists.r-forge.r-project.org/pipermail/rcpp-devel/2021-July/010648.html)
#' thread appears to provide some insight. In a fresh session, run
#' `remove.packages(c("MatchIt", "Rcpp")); install.packages("MatchIt")`.
#' This should sync *MatchIt* and *Rcpp* and ensure they work
#' correctly.
#'
#' @section Outputs:
#' All outputs described in [matchit()] are returned with
#' `method = "nearest"`. When `replace = TRUE`, the `subclass`
#' component is omitted. `include.obj` is ignored.
#'
#' @details
#' ## Mahalanobis Distance Matching
#'
#' Mahalanobis distance matching can be done one of two ways:
#' \enumerate{
#' \item{If no propensity score needs to be estimated, `distance` should be
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
#' \item{If a propensity score needs to be estimated for any reason, e.g., for
#' common support with `discard` or for creating a caliper,
#' `distance` should be whatever method is used to estimate the propensity
#' score or a vector of distance measures. Use `mahvars` to specify the
#' variables used to create the Mahalanobis distance. For example, to perform
#' Mahalanobis within a propensity score caliper, the following could be run:
#'
#' \preformatted{
#' matchit(treat ~ X1 + X2 + X3, method = "nearest",
#'         distance = "glm", caliper = .25,
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
#' ## Estimand
#'
#' The `estimand` argument controls whether control units are selected to be
#' matched with treated units (`estimand = "ATT"`) or treated units are
#' selected to be matched with control units (`estimand = "ATC"`). The
#' "focal" group (e.g., the treated units for the ATT) is typically made to be
#' the smaller treatment group, and a warning will be thrown if it is not set
#' that way unless `replace = TRUE`. Setting `estimand = "ATC"` is
#' equivalent to swapping all treated and control labels for the treatment
#' variable. When `estimand = "ATC"`, the default `m.order` is
#' `"smallest"`, and the `match.matrix` component of the output will
#' have the names of the control units as the rownames and be filled with the
#' names of the matched treated units (opposite to when `estimand =
#' "ATT"`). Note that the argument supplied to `estimand` doesn't
#' necessarily correspond to the estimand actually targeted; it is merely a
#' switch to trigger which treatment group is considered "focal".
#'
#' ## Variable Ratio Matching
#'
#' `matchit()` can perform variable
#' ratio "extremal" matching as described by Ming and Rosenbaum (2000). This
#' method tends to result in better balance than fixed ratio matching at the
#' expense of some precision. When `ratio > 1`, rather than requiring all
#' treated units to receive `ratio` matches, each treated unit is assigned
#' a value that corresponds to the number of control units they will be matched
#' to. These values are controlled by the arguments `min.controls` and
#' `max.controls`, which correspond to \eqn{\alpha} and \eqn{\beta},
#' respectively, in Ming and Rosenbaum (2000), and trigger variable ratio
#' matching to occur. Some treated units will receive `min.controls`
#' matches and others will receive `max.controls` matches (and one unit
#' may have an intermediate number of matches); how many units are assigned
#' each number of matches is determined by the algorithm described in Ming and
#' Rosenbaum (2000, p119). `ratio` controls how many total control units
#' will be matched: `n1 * ratio` control units will be matched, where
#' `n1` is the number of treated units, yielding the same total number of
#' matched controls as fixed ratio matching does.
#'
#' Variable ratio matching cannot be used with Mahalanobis distance matching or
#' when `distance` is supplied as a matrix. The calculations of the
#' numbers of control units each treated unit will be matched to occurs without
#' consideration of `caliper` or `discard`. `ratio` does not
#' have to be an integer but must be greater than 1 and less than `n0/n1`,
#' where `n0` and `n1` are the number of control and treated units,
#' respectively. Setting `ratio = n0/n1` performs a crude form of full
#' matching where all control units are matched. If `min.controls` is not
#' specified, it is set to 1 by default. `min.controls` must be less than
#' `ratio`, and `max.controls` must be greater than `ratio`. See
#' Examples below for an example of their use.
#'
#' @seealso [matchit()] for a detailed explanation of the inputs and outputs of
#' a call to `matchit()`.
#'
#' [method_optimal()] for optimal pair matching, which is similar to
#' nearest neighbor matching except that an overall distance criterion is
#' minimized.
#'
#' @references In a manuscript, you don't need to cite another package when
#' using `method = "nearest"` because the matching is performed completely
#' within *MatchIt*. For example, a sentence might read:
#'
#' *Nearest neighbor matching was performed using the MatchIt package (Ho, Imai, King, & Stuart, 2011) in R.*
#'
#' @examples
#' data("lalonde")
#'
#' # 1:1 greedy NN matching on the PS
#' m.out1 <- matchit(treat ~ age + educ + race + nodegree +
#'                     married + re74 + re75, data = lalonde,
#'                   method = "nearest")
#' m.out1
#' summary(m.out1)
#'
#' # 3:1 NN Mahalanobis distance matching with
#' # replacement within a PS caliper
#' m.out2 <- matchit(treat ~ age + educ + race + nodegree +
#'                     married + re74 + re75, data = lalonde,
#'                   method = "nearest", replace = TRUE,
#'                   mahvars = ~ age + educ + re74 + re75,
#'                   ratio = 3, caliper = .02)
#' m.out2
#' summary(m.out2, un = FALSE)
#'
#' # 1:1 NN Mahalanobis distance matching within calipers
#' # on re74 and re75 and exact matching on married and race
#' m.out3 <- matchit(treat ~ age + educ + re74 + re75, data = lalonde,
#'                   method = "nearest", distance = "mahalanobis",
#'                   exact = ~ married + race,
#'                   caliper = c(re74 = .2, re75 = .15))
#' m.out3
#' summary(m.out3, un = FALSE)
#'
#' # 2:1 variable ratio NN matching on the PS
#' m.out4 <- matchit(treat ~ age + educ + race + nodegree +
#'                     married + re74 + re75, data = lalonde,
#'                   method = "nearest", ratio = 2,
#'                   min.controls = 1, max.controls = 12)
#' m.out4
#' summary(m.out4, un = FALSE)
#'
#' # Some units received 1 match and some received 12
#' table(table(m.out4$subclass[m.out4$treat == 0]))
#'
NULL

matchit2nearest <- function(treat, data, distance, discarded,
                             ratio = 1, s.weights = NULL, replace = FALSE, m.order = NULL,
                             caliper = NULL, mahvars = NULL, exact = NULL,
                             formula = NULL, estimand = "ATT", verbose = FALSE,
                             is.full.mahalanobis,
                             antiexact = NULL, unit.id = NULL, ...){

  if (verbose) {
    rlang::check_installed("RcppProgress")
    cat("Nearest neighbor matching... \n")
  }

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

  treat <- setNames(as.integer(treat == focal), names(treat))

  if (is.full.mahalanobis) {
    if (length(attr(terms(formula, data = data), "term.labels")) == 0) {
      .err(sprintf("covariates must be specified in the input formula when `distance = \"%s\"`",
                   attr(is.full.mahalanobis, "transform")))
    }
    mahvars <- formula
  }

  n.obs <- length(treat)
  n1 <- sum(treat == 1)
  n0 <- n.obs - n1

  lab <- names(treat)
  lab1 <- lab[treat == 1]

  if (!is.null(distance)) {
    names(distance) <- names(treat)
  }

  min.controls <- attr(ratio, "min.controls")
  max.controls <- attr(ratio, "max.controls")

  mahcovs <- distance_mat <- NULL
  if (!is.null(mahvars)) {
    transform <- if (is.full.mahalanobis) attr(is.full.mahalanobis, "transform") else "mahalanobis"
    mahcovs <- transform_covariates(mahvars, data = data, method = transform,
                                    s.weights = s.weights, treat = treat,
                                    discarded = discarded)
  }
  else if (is.matrix(distance)) {
    distance_mat <- distance
    distance <- NULL
  }

  #Process caliper
  if (!is.null(caliper)) {
    if (any(names(caliper) != "")) {
      caliper.covs <- caliper[names(caliper) != ""]
      caliper.covs.mat <- get.covs.matrix(reformulate(names(caliper.covs)), data = data)
    }
    else {
      caliper.covs.mat <- caliper.covs <- NULL
    }

    if (any(names(caliper) == "")) {
      caliper.dist <- caliper[names(caliper) == ""]
    }
    else {
      caliper.dist <- NULL
    }
  }
  else {
    caliper.dist <- caliper.covs <- NULL
    caliper.covs.mat <- NULL
  }

  #Process antiexact
  if (!is.null(antiexact)) {
    antiexactcovs <- model.frame(antiexact, data)
    antiexactcovs <- do.call("cbind", lapply(seq_len(ncol(antiexactcovs)), function(i) {
      as.integer(as.factor(antiexactcovs[[i]]))
    }))
  }
  else {
    antiexactcovs <- NULL
  }

  reuse.max <- attr(replace, "reuse.max")

  if (reuse.max >= n1) {
    m.order <- "data"
  }

  if (!is.null(unit.id) && reuse.max < n1) {
    unit.id <- process.variable.input(unit.id, data)
    unit.id <- factor(exactify(model.frame(unit.id, data = data),
                               nam = lab, sep = ", ", include_vars = TRUE))
    num_ctrl_unit.ids <- length(unique(unit.id[treat == 0]))

    #If each control unit is a unit.id, unit.ids are meaningless
    if (num_ctrl_unit.ids == n0) unit.id <- NULL
  }
  else {
    unit.id <- NULL
  }

  if (!is.null(exact)) {
    ex <- factor(exactify(model.frame(exact, data = data), nam = lab, sep = ", ", include_vars = TRUE))

    cc <- intersect(as.integer(ex)[treat==1], as.integer(ex)[treat==0])
    if (length(cc) == 0) .err("No matches were found")

    if (reuse.max < n1) {

      e_ratios <- vapply(levels(ex), function(e) {
        if (is.null(unit.id)) sum(treat[ex == e] == 0)*(reuse.max/sum(treat[ex == e] == 1))
        else length(unique(unit.id[treat == 0 & ex == e]))*(reuse.max/sum(treat[ex == e] == 1))
      }, numeric(1L))

      if (any(e_ratios < 1)) {
        .wrn(sprintf("fewer %s units than %s units in some `exact` strata; not all %s units will get a match",
                        tc[2], tc[1], tc[1]))
      }
      if (ratio > 1 && any(e_ratios < ratio)) {
        if (is.null(max.controls) || ratio == max.controls)
          .wrn(sprintf("not all %s units will get %s matches",
                          tc[1], ratio))
        else
          .wrn(sprintf("not enough %s units for an average of %s matches per %s unit in all `exact` strata",
                          tc[2], ratio, tc[1]))
      }
    }
  }
  else {
    ex <- NULL

    if (reuse.max < n1) {

      e_ratios <- {
        if (is.null(unit.id)) as.numeric(reuse.max)*n0/n1
        else as.numeric(reuse.max)*num_ctrl_unit.ids/n1
      }

      if (e_ratios < 1) {
        .wrn(sprintf("fewer %s %s than %s units; not all %s units will get a match",
                        tc[2], if (is.null(unit.id)) "units" else "unit IDs", tc[1], tc[1]))
      }
      else if (e_ratios < ratio) {
        if (is.null(max.controls) || ratio == max.controls)
          .wrn(sprintf("not all %s units will get %s matches",
                          tc[1], ratio))
        else
          .wrn(sprintf("not enough %s %s for an average of %s matches per %s unit",
                          tc[2], if (is.null(unit.id)) "units" else "unit IDs", ratio, tc[1]))
      }
    }
  }

  #Variable ratio (extremal matching), Ming & Rosenbaum (2000)
  #Each treated unit get its own value of ratio
  if (!is.null(max.controls)) {
    if (is.null(distance)) {
      if (is.full.mahalanobis) .err(sprintf("`distance` cannot be \"%s\" for variable ratio matching",
                                            transform))
      .err("`distance` cannot be supplied as a matrix for variable ratio matching")
    }

    m <- round(ratio * n1)
    # if (m > sum(treat == 0)) stop("'ratio' must be less than or equal to n0/n1")

    kmax <- floor((m - min.controls*(n1-1)) / (max.controls - min.controls))
    kmin <- n1 - kmax - 1
    kmed <- m - (min.controls*kmin + max.controls*kmax)

    ratio0 <- c(rep(min.controls, kmin), kmed, rep(max.controls, kmax))

    #Make sure no units are assigned 0 matches
    if (any(ratio0 == 0)) {
      ind <- which(ratio0 == 0)
      ratio0[ind] <- 1
      ratio0[ind + 1] <- ratio0[ind + 1] - 1
    }

    ratio <- rep(NA_integer_, n1)

    #Order by distance; treated are supposed to have higher values
    ratio[order(distance[treat == 1],
                decreasing = mean(distance[treat == 1]) > mean(distance[treat != 1]))] <- ratio0
    ratio <- as.integer(ratio)
  }
  else {
    ratio <- as.integer(rep(ratio, n1))
  }

  m.order <- {
    if (is.null(distance)) match_arg(m.order, c("data", "random"))
    else if (!is.null(m.order)) match_arg(m.order, c("largest", "smallest", "random", "data"))
    else if (estimand == "ATC") "smallest"
    else "largest"
  }

  if (is.null(ex) || !is.null(unit.id)) {
    ord <- switch(m.order,
                  "largest" = order(distance[treat == 1], decreasing = TRUE),
                  "smallest" = order(distance[treat == 1], decreasing = FALSE),
                  "random" = sample.int(n1),
                  "data" = seq_len(n1))

    mm <- nn_matchC_dispatch(treat, ord, ratio, discarded, reuse.max, distance, distance_mat,
                             ex, caliper.dist, caliper.covs, caliper.covs.mat, mahcovs,
                             antiexactcovs, unit.id, verbose)

  }
  else {
    distance_ <- caliper.covs.mat_ <- mahcovs_ <- antiexactcovs_ <- distance_mat_ <- NULL
    mm_list <- lapply(levels(ex)[cc], function(e) {
      if (verbose) {
        cat(sprintf("Matching subgroup %s/%s: %s...\n",
                    match(e, levels(ex)), nlevels(ex), e))
      }

      .e <- which(ex == e)
      treat_ <- treat[.e]
      discarded_ <- discarded[.e]
      if (!is.null(distance)) distance_ <- distance[.e]
      if (!is.null(caliper.covs.mat)) caliper.covs.mat_ <- caliper.covs.mat[.e,,drop = FALSE]
      if (!is.null(mahcovs)) mahcovs_ <- mahcovs[.e,,drop = FALSE]
      if (!is.null(antiexactcovs)) antiexactcovs_ <- antiexactcovs[.e,,drop = FALSE]
      if (!is.null(distance_mat)) {
        .e1 <- which(ex[treat==1] == e)
        .e0 <- which(ex[treat==0] == e)
        distance_mat_ <- distance_mat[.e1, .e0, drop = FALSE]
      }
      ratio_ <- ratio[ex[treat==1]==e]

      n1_ <- sum(treat_ == 1)
      ord_ <- switch(m.order,
                     "largest" = order(distance_[treat_ == 1], decreasing = TRUE),
                     "smallest" = order(distance_[treat_ == 1], decreasing = FALSE),
                     "random" = sample.int(n1_),
                     "data" = seq_len(n1_))

      mm_ <- nn_matchC_dispatch(treat_, ord_, ratio_, discarded_, reuse.max, distance_, distance_mat_,
                               NULL, caliper.dist, caliper.covs, caliper.covs.mat_, mahcovs_,
                               antiexactcovs_, NULL, verbose)

      #Ensure matched indices correspond to indices in full sample, not subgroup
      mm_[] <- seq_along(treat)[.e][mm_]
      mm_
    })

    #Construct match.matrix
    mm <- matrix(NA_integer_, nrow = length(lab1),
                 ncol = max(vapply(mm_list, ncol, numeric(1L))),
                 dimnames = list(lab1, NULL))

    for (m in mm_list) {
      mm[rownames(m), seq_len(ncol(m))] <- m
    }
  }

  if (verbose) cat("Calculating matching weights... ")

  if (reuse.max > 1) {
    psclass <- NULL
  }
  else {
    psclass <- mm2subclass(mm, treat)
  }

  res <- list(match.matrix = nummm2charmm(mm, treat),
              subclass = psclass,
              weights = weights.matrix(mm, treat))

  if (verbose) cat("Done.\n")

  class(res) <- "matchit"

  res
}

# Dispatches Rcpp function for NN matching
# nn_matchC_vec() if distance_mat and mahcovs are NULL
# nn_matchC() otherwise
nn_matchC_dispatch <- function(treat, ord, ratio, discarded, reuse.max, distance, distance_mat, ex, caliper.dist,
                               caliper.covs, caliper.covs.mat, mahcovs, antiexactcovs, unit.id, verbose) {
  if (is.null(distance_mat) && is.null(mahcovs)) {
    nn_matchC_vec(treat, ord, ratio, discarded, reuse.max, distance,
                  ex, caliper.dist, caliper.covs, caliper.covs.mat,
                  antiexactcovs, unit.id, verbose)
  }
  else {
    nn_matchC(treat, ord, ratio, discarded, reuse.max, distance, distance_mat,
              ex, caliper.dist, caliper.covs, caliper.covs.mat, mahcovs,
              antiexactcovs, unit.id, verbose)
  }
}