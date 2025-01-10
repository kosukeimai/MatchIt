#' Genetic Matching
#' @name method_genetic
#' @aliases method_genetic
#' @usage NULL
#'
#' @description
#' In [matchit()], setting `method = "genetic"` performs genetic matching.
#' Genetic matching is a form of nearest neighbor matching where distances are
#' computed as the generalized Mahalanobis distance, which is a generalization
#' of the Mahalanobis distance with a scaling factor for each covariate that
#' represents the importance of that covariate to the distance. A genetic
#' algorithm is used to select the scaling factors. The scaling factors are
#' chosen as those which maximize a criterion related to covariate balance,
#' which can be chosen, but which by default is the smallest p-value in
#' covariate balance tests among the covariates. This method relies on and is a
#' wrapper for \pkgfun{Matching}{GenMatch} and \pkgfun{Matching}{Match}, which use
#' \pkgfun{rgenoud}{genoud} to perform the optimization using the genetic
#' algorithm.
#'
#' This page details the allowable arguments with `method = "genetic"`.
#' See [matchit()] for an explanation of what each argument means in a general
#' context and how it can be specified.
#'
#' Below is how `matchit()` is used for genetic matching:
#' \preformatted{
#' matchit(formula,
#'         data = NULL,
#'         method = "genetic",
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
#'         replace = FALSE,
#'         m.order = NULL,
#'         caliper = NULL,
#'         ratio = 1,
#'         verbose = FALSE,
#'         ...) }
#'
#' @param formula a two-sided [formula] object containing the treatment and
#' covariates to be used in creating the distance measure used in the matching.
#' This formula will be supplied to the functions that estimate the distance
#' measure and is used to determine the covariates whose balance is to be
#' optimized.
#' @param data a data frame containing the variables named in `formula`.
#' If not found in `data`, the variables will be sought in the
#' environment.
#' @param method set here to `"genetic"`.
#' @param distance the distance measure to be used. See [`distance`]
#' for allowable options. When set to a method of estimating propensity scores
#' or a numeric vector of distance values, the distance measure is included
#' with the covariates in `formula` to be supplied to the generalized
#' Mahalanobis distance matrix unless `mahvars` is specified. Otherwise,
#' only the covariates in `formula` are supplied to the generalized
#' Mahalanobis distance matrix to have their scaling factors chosen.
#' `distance` *cannot* be supplied as a distance matrix. Supplying
#' any method of computing a distance matrix (e.g., `"mahalanobis"`) has
#' the same effect of omitting propensity score but does not affect how the
#' distance between units is computed otherwise.
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
#' @param mahvars when a distance corresponds to a propensity score (e.g., for
#' caliper matching or to discard units for common support), which covariates
#' should be supplied to the generalized Mahalanobis distance matrix for
#' matching. If unspecified, all variables in `formula` will be supplied
#' to the distance matrix. Use `mahvars` to only supply a subset. Even if
#' `mahvars` is specified, balance will be optimized on all covariates in
#' `formula`. See Details.
#' @param antiexact for which variables anti-exact matching should take place.
#' Anti-exact matching is processed using the `restrict` argument to
#' `Matching::GenMatch()` and `Matching::Match()`.
#' @param discard a string containing a method for discarding units outside a
#' region of common support. Only allowed when `distance` corresponds to a
#' propensity score.
#' @param reestimate if `discard` is not `"none"`, whether to
#' re-estimate the propensity score in the remaining sample prior to matching.
#' @param s.weights the variable containing sampling weights to be incorporated
#' into propensity score models and balance statistics. These are also supplied
#' to `GenMatch()` for use in computing the balance t-test p-values in the
#' process of matching.
#' @param replace whether matching should be done with replacement.
#' @param m.order the order that the matching takes place. Allowable options
#'   include `"largest"`, where matching takes place in descending order of
#'   distance measures; `"smallest"`, where matching takes place in ascending
#'   order of distance measures; `"random"`, where matching takes place
#'   in a random order; and `"data"` where matching takes place based on the
#'   order of units in the data. When `m.order = "random"`, results may differ
#'   across different runs of the same code unless a seed is set and specified
#'   with [set.seed()]. The default of `NULL` corresponds to `"largest"` when a
#'   propensity score is estimated or supplied as a vector and `"data"`
#'   otherwise.
#' @param caliper the width(s) of the caliper(s) used for caliper matching. See
#' Details and Examples.
#' @param std.caliper `logical`; when calipers are specified, whether they
#' are in standard deviation units (`TRUE`) or raw units (`FALSE`).
#' @param ratio how many control units should be matched to each treated unit
#' for k:1 matching. Should be a single integer value.
#' @param verbose `logical`; whether information about the matching
#' process should be printed to the console. When `TRUE`, output from
#' `GenMatch()` with `print.level = 2` will be displayed. Default is
#' `FALSE` for no printing other than warnings.
#' @param \dots additional arguments passed to \pkgfun{Matching}{GenMatch}.
#' Potentially useful options include `pop.size`, `max.generations`,
#' and `fit.func`. If `pop.size` is not specified, a warning from
#' *Matching* will be thrown reminding you to change it. Note that the
#' `ties` and `CommonSupport` arguments are set to `FALSE` and
#' cannot be changed. If `distance.tolerance` is not specified, it is set
#' to 0, whereas the default in *Matching* is 1e-5.
#'
#' @section Outputs:
#' All outputs described in [matchit()] are returned with
#' `method = "genetic"`. When `replace = TRUE`, the `subclass`
#' component is omitted. When `include.obj = TRUE` in the call to
#' `matchit()`, the output of the call to \pkgfun{Matching}{GenMatch} will be
#' included in the output.
#'
#' @details
#' In genetic matching, covariates play three roles: 1) as the variables on
#' which balance is optimized, 2) as the variables in the generalized
#' Mahalanobis distance between units, and 3) in estimating the propensity
#' score. Variables supplied to `formula` are always used for role (1), as
#' the variables on which balance is optimized. When `distance`
#' corresponds to a propensity score, the covariates are also used to estimate
#' the propensity score (unless it is supplied). When `mahvars` is
#' specified, the named variables will form the covariates that go into the
#' distance matrix. Otherwise, the variables in `formula` along with the
#' propensity score will go into the distance matrix. This leads to three ways
#' to use `distance` and `mahvars` to perform the matching:
#'
#' \enumerate{
#' \item{When `distance` corresponds to a propensity score and `mahvars`
#' *is not* specified, the covariates in `formula` along with the
#' propensity score are used to form the generalized Mahalanobis distance
#' matrix. This is the default and most typical use of `method =
#' "genetic"` in `matchit()`.
#' }
#' \item{When `distance` corresponds to a propensity score and `mahvars`
#' *is* specified, the covariates in `mahvars` are used to form the
#' generalized Mahalanobis distance matrix. The covariates in `formula`
#' are used to estimate the propensity score and have their balance optimized
#' by the genetic algorithm. The propensity score is not included in the
#' generalized Mahalanobis distance matrix.
#' }
#' \item{When `distance` is a method of computing a distance matrix
#' (e.g.,`"mahalanobis"`), no propensity score is estimated, and the
#' covariates in `formula` are used to form the generalized Mahalanobis
#' distance matrix. Which specific method is supplied has no bearing on how the
#' distance matrix is computed; it simply serves as a signal to omit estimation
#' of a propensity score.
#' }
#' }
#'
#' When a caliper is specified, any variables mentioned in `caliper`,
#' possibly including the propensity score, will be added to the matching
#' variables used to form the generalized Mahalanobis distance matrix. This is
#' because *Matching* doesn't allow for the separation of caliper
#' variables and matching variables in genetic matching.
#'
#' ## Estimand
#'
#' The `estimand` argument controls whether control
#' units are selected to be matched with treated units (`estimand =
#' "ATT"`) or treated units are selected to be matched with control units
#' (`estimand = "ATC"`). The "focal" group (e.g., the treated units for
#' the ATT) is typically made to be the smaller treatment group, and a warning
#' will be thrown if it is not set that way unless `replace = TRUE`.
#' Setting `estimand = "ATC"` is equivalent to swapping all treated and
#' control labels for the treatment variable. When `estimand = "ATC"`, the
#' default `m.order` is `"smallest"`, and the `match.matrix`
#' component of the output will have the names of the control units as the
#' rownames and be filled with the names of the matched treated units (opposite
#' to when `estimand = "ATT"`). Note that the argument supplied to
#' `estimand` doesn't necessarily correspond to the estimand actually
#' targeted; it is merely a switch to trigger which treatment group is
#' considered "focal". Note that while `GenMatch()` and `Match()`
#' support the ATE as an estimand, `matchit()` only supports the ATT and
#' ATC for genetic matching.
#'
#' ## Reproducibility
#'
#' Genetic matching involves a random component, so a seed must be set using [set.seed()] to ensure reproducibility. When `cluster` is used for parallel processing, the seed must be compatible with parallel processing (e.g., by setting `type = "L'Ecuyer-CMRG"`).
#'
#' @seealso [matchit()] for a detailed explanation of the inputs and outputs of
#' a call to `matchit()`.
#'
#' \pkgfun{Matching}{GenMatch} and \pkgfun{Matching}{Match}, which do the work.
#'
#' @references In a manuscript, be sure to cite the following papers if using
#' `matchit()` with `method = "genetic"`:
#'
#' Diamond, A., & Sekhon, J. S. (2013). Genetic matching for estimating causal
#' effects: A general multivariate matching method for achieving balance in
#' observational studies. Review of Economics and Statistics, 95(3), 932–945. \doi{10.1162/REST_a_00318}
#'
#' Sekhon, J. S. (2011). Multivariate and Propensity Score Matching Software
#' with Automated Balance Optimization: The Matching package for R. Journal of
#' Statistical Software, 42(1), 1–52. \doi{10.18637/jss.v042.i07}
#'
#' For example, a sentence might read:
#'
#' *Genetic matching was performed using the MatchIt package (Ho, Imai,
#' King, & Stuart, 2011) in R, which calls functions from the Matching package
#' (Diamond & Sekhon, 2013; Sekhon, 2011).*
#'
#' @examplesIf all(sapply(c("Matching", "rgenoud"), requireNamespace, quietly = TRUE))
#' data("lalonde")
#'
#' # 1:1 genetic matching with PS as a covariate
#' m.out1 <- matchit(treat ~ age + educ + race + nodegree +
#'                     married + re74 + re75, data = lalonde,
#'                   method = "genetic",
#'                   pop.size = 10) #use much larger pop.size
#' m.out1
#' summary(m.out1)
#'
#' # 2:1 genetic matching with replacement without PS
#' m.out2 <- matchit(treat ~ age + educ + race + nodegree +
#'                     married + re74 + re75, data = lalonde,
#'                   method = "genetic", replace = TRUE,
#'                   ratio = 2, distance = "mahalanobis",
#'                   pop.size = 10) #use much larger pop.size
#' m.out2
#' summary(m.out2, un = FALSE)
#'
#' # 1:1 genetic matching on just age, educ, re74, and re75
#' # within calipers on PS and educ; other variables are
#' # used to estimate PS
#' m.out3 <- matchit(treat ~ age + educ + race + nodegree +
#'                     married + re74 + re75, data = lalonde,
#'                   method = "genetic",
#'                   mahvars = ~ age + educ + re74 + re75,
#'                   caliper = c(.05, educ = 2),
#'                   std.caliper = c(TRUE, FALSE),
#'                   pop.size = 10) #use much larger pop.size
#' m.out3
#' summary(m.out3, un = FALSE)
NULL

matchit2genetic <- function(treat, data, distance, discarded,
                            ratio = 1, s.weights = NULL, replace = FALSE, m.order = NULL,
                            caliper = NULL, mahvars = NULL, exact = NULL,
                            formula = NULL, estimand = "ATT", verbose = FALSE,
                            is.full.mahalanobis, use.genetic = TRUE,
                            antiexact = NULL, ...) {

  rlang::check_installed(c("Matching", "rgenoud"))

  .cat_verbose("Genetic matching...\n", verbose = verbose)

  args <- names(formals(Matching::GenMatch))
  A <- ...mget(args)
  A[lengths(A) == 0L] <- NULL

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

  if (!replace) {
    if (sum(!discarded & treat != focal) < sum(!discarded & treat == focal)) {
      .wrn(sprintf("fewer %s units than %s units; not all %s units will get a match",
                   tc[2], tc[1], tc[1]))
    }
    else if (sum(!discarded & treat != focal) < sum(!discarded & treat == focal)*ratio) {
      .err(sprintf("not enough %s units for %s matches for each %s unit",
                   tc[2], ratio, tc[1]))
    }
  }

  treat <- setNames(as.integer(treat == focal), names(treat))

  n.obs <- length(treat)
  n1 <- sum(treat == 1)

  if (is_null(names(treat))) names(treat) <- seq_len(n.obs)

  m.order <- {
    if (is_null(distance)) match_arg(m.order, c("data", "random"))
    else if (is_not_null(m.order)) match_arg(m.order, c("largest", "smallest", "random", "data"))
    else if (estimand == "ATC") "smallest"
    else "largest"
  }

  ord <- switch(m.order,
                "largest" = order(distance, decreasing = TRUE),
                "smallest" = order(distance),
                "random" = sample.int(n.obs),
                "data" = seq_len(n.obs))
  ord <- ord[!ord %in% which(discarded)]

  #Create X (matching variables) and covs_to_balance
  covs_to_balance <- get_covs_matrix(formula, data = data)
  if (is_not_null(mahvars)) {
    X <- get_covs_matrix_for_dist(mahvars, data = data)
  }
  else if (is.full.mahalanobis) {
    X <- covs_to_balance
  }
  else {
    X <- cbind(covs_to_balance, distance)
  }

  if (ncol(covs_to_balance) == 0L) {
    .err("covariates must be specified in the input formula to use genetic matching")
  }

  #Process exact; exact.log will be supplied to GenMatch() and Match()
  if (is_not_null(exact)) {
    #Add covariates in exact not in X to X
    ex <- unclass(exactify(model.frame(exact, data = data), names(treat),
                           sep = ", ", include_vars = TRUE))

    cc <- intersect(ex[treat==1], ex[treat==0])

    if (is_null(cc)) {
      .err("No matches were found")
    }

    X <- cbind(X, ex)

    exact.log <- c(rep.int(FALSE, ncol(X) - 1L), TRUE)
  }
  else {
    exact.log <- ex <- NULL
  }

  #Reorder data according to m.order since Match matches in order of data;
  #ord already excludes discarded units

  treat_ <- treat[ord]
  covs_to_balance <- covs_to_balance[ord,,drop = FALSE]
  X <- X[ord,,drop = FALSE]
  if (is_not_null(s.weights)) s.weights <- s.weights[ord]

  #Process caliper; cal will be supplied to GenMatch() and Match()
  cal <- dist.cal <- cov.cals <- NULL
  if (is_not_null(caliper) && any(caliper < 0)) {
    neg.cal <- names(caliper)[caliper < 0]

    if (any(neg.cal != "")) {
      negcalcovs <- get_covs_matrix(reformulate(neg.cal[neg.cal != ""]), data = data)[ord,,drop = FALSE]
      negcalcovs_restrict <- do.call("rbind", lapply(seq_len(ncol(negcalcovs)), function(i) {
        do.call("rbind", lapply(which(treat_ == 1), function(j) {
          restricted_controls <- which(treat_ == 0 & abs(negcalcovs[j, i] - negcalcovs[, i]) <= -caliper[neg.cal[neg.cal != ""][i]])

          if (is_null(restricted_controls)) {
            return(NULL)
          }

          cbind(j, restricted_controls, -1)
        }))
      }))

      if (is_not_null(negcalcovs_restrict)) {
        A[["restrict"]] <- {
          if (is_null(A[["restrict"]])) unique(negcalcovs_restrict)
          else rbind(A[["restrict"]], unique(negcalcovs_restrict))
        }
      }
    }

    if (any(neg.cal == "")) {
      negcaldist_restrict <- do.call("rbind", lapply(which(treat_ == 1), function(j) {
        restricted_controls <- which(treat_ == 0 & abs(distance[ord][j] - distance[ord]) <= -caliper[names(caliper) == ""])

        if (is_null(restricted_controls)) {
          return(NULL)
        }

        cbind(j, restricted_controls, -1)
      }))

      if (is_not_null(negcaldist_restrict)) {
        A[["restrict"]] <- {
          if (is_null(A[["restrict"]])) unique(negcaldist_restrict)
          else rbind(A[["restrict"]], unique(negcaldist_restrict))
        }
      }
    }

    caliper <- caliper[caliper >= 0]
  }

  #Add covariates in caliper other than distance (cov.cals) not in X to X
  if (is_not_null(caliper)) {
    cov.cals <- setdiff(names(caliper), "")
    if (is_not_null(cov.cals) && !all(cov.cals %in% colnames(X))) {
      calcovs <- get_covs_matrix(reformulate(cov.cals[!cov.cals %in% colnames(X)]), data = data)[ord,,drop = FALSE]
      X <- cbind(X, calcovs)

      #Expand exact.log for newly added covariates
      if (is_not_null(exact.log)) {
        exact.log <- c(exact.log, rep.int(FALSE, ncol(calcovs)))
      }
    }

    #Matching::Match multiplies calipers by pop SD, so we need to divide by pop SD to unstandardize
    pop.sd <- function(x) sqrt(sum((x-mean(x))^2)/length(x))
    caliper <- caliper / vapply(names(caliper), function(x) {
      if (x == "") pop.sd(distance[ord])
      else pop.sd(X[, x])
    }, numeric(1L))

    #cal needs one value per variable in X
    cal <- setNames(rep.int(Inf, ncol(X)), colnames(X))

    #First put covariate calipers into cal
    if (is_not_null(cov.cals)) {
      cal[intersect(cov.cals, names(cal))] <- caliper[intersect(cov.cals, names(cal))]
    }

    #Then put distance caliper into cal
    if (hasName(caliper, "")) {
      dist.cal <- caliper[names(caliper) == ""]
      if (is_not_null(mahvars)) {
        #If mahvars specified, distance is not yet in X, so add it to X
        X <- cbind(X, distance[ord])
        cal <- c(cal, dist.cal)

        #Expand exact.log for newly added distance
        if (is_not_null(exact.log)) exact.log <- c(exact.log, FALSE)
      }
      else {
        #Otherwise, distance is in X at the specified index
        cal[ncol(covs_to_balance) + 1] <- dist.cal
      }
    }
    else {
      dist.cal <- NULL
    }
  }

  if (is_not_null(antiexact)) {
    antiexactcovs <- model.frame(antiexact, data)[ord,,drop = FALSE]
    antiexact_restrict <- do.call("rbind", lapply(seq_len(ncol(antiexactcovs)), function(i) {
      do.call("rbind", lapply(which(treat_ == 1), function(j) {
        restricted_controls <- which(treat_ == 0 & antiexactcovs[[i]][j] == antiexactcovs[[i]])

        if (is_null(restricted_controls)) {
          return(NULL)
        }

        cbind(j, restricted_controls, -1)
      }))
    }))

    if (is_not_null(antiexact_restrict)) {
      A[["restrict"]] <- {
        if (is_null(A[["restrict"]])) unique(antiexact_restrict)
        else rbind(A[["restrict"]], unique(antiexact_restrict))
      }
    }
  }
  else {
    antiexactcovs <- NULL
  }

  if (is_null(A[["distance.tolerance"]])) {
    A[["distance.tolerance"]] <- 0
  }

  if (use.genetic) {
      matchit_try({
        g.out <- do.call(Matching::GenMatch,
                         c(list(Tr = treat_, X = X, BalanceMatrix = covs_to_balance,
                                M = ratio, exact = exact.log, caliper = cal,
                                replace = replace, estimand = "ATT", ties = FALSE,
                                CommonSupport = FALSE, verbose = verbose,
                                weights = s.weights, print.level = 2*verbose),
                           A[names(A) %in% args]))
      }, from = "Matching",
      dont_warn_if = c("replace==FALSE, but there are more (weighted) treated obs than control obs",
                       "no valid matches"))
  }
  else {
    #For debugging
    g.out <- NULL
  }

  lab <- names(treat)
  lab1 <- lab[treat == 1]

  lab_ <- names(treat_)

  ind_ <- seq_along(treat)[ord]

  matchit_try({
    m.out <- Matching::Match(Tr = treat_, X = X,
                             M = ratio, exact = exact.log, caliper = cal,
                             replace = replace, estimand = "ATT", ties = FALSE,
                             weights = s.weights, CommonSupport = FALSE,
                             distance.tolerance = A[["distance.tolerance"]], Weight = 3,
                             Weight.matrix = {
                               if (use.genetic) g.out
                               else if (is_null(s.weights)) generalized_inverse(cor(X))
                               else generalized_inverse(cov.wt(X, s.weights, cor = TRUE)$cor)
                             },
                             restrict = A[["restrict"]], version = "fast")
  }, from = "Matching",
  dont_warn_if = c("replace==FALSE, but there are more (weighted) treated obs than control obs",
                   "no valid matches"))

  if (typeof(m.out) == "logical" && all(is.na(m.out))) {
    .err("no units were matched")
  }

  #Note: must use character match.matrix because of re-ordering treat into treat_
  mm <- matrix(NA_integer_, nrow = n1, ncol = max(table(m.out$index.treated)),
               dimnames = list(lab1, NULL))

  unique.matched.focal <- unique(m.out$index.treated, nmax = n1)

  ind1__ <- match(lab_, lab1)
  for (i in unique.matched.focal) {
    matched.units <- ind_[m.out$index.control[m.out$index.treated == i]]
    mm[ind1__[i], seq_along(matched.units)] <- matched.units
  }

  .cat_verbose("Calculating matching weights... ", verbose = verbose)

  if (replace) {
    psclass <- NULL
    weights <- get_weights_from_mm(mm, treat, 1L)
  }
  else {
    psclass <- mm2subclass(mm, treat, 1L)
    weights <- get_weights_from_subclass(psclass, treat)
  }

  res <- list(match.matrix = nummm2charmm(mm, treat),
              subclass = psclass,
              weights = weights,
              obj = g.out)

  .cat_verbose("Done.\n", verbose = verbose)

  class(res) <- "matchit"
  res
}
