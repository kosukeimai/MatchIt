#' View a balance summary of a `matchit` object
#'
#' Computes and prints balance statistics for `matchit` and
#' `matchit.subclass` objects. Balance should be assessed to ensure the
#' matching or subclassification was effective at eliminating treatment group
#' imbalance and should be reported in the write-up of the results of the
#' analysis.
#'
#' @aliases summary.matchit summary.matchit.subclass print.summary.matchit
#' print.summary.matchit.subclass
#'
#' @param object a `matchit` object; the output of a call to [matchit()].
#' @param interactions `logical`; whether to compute balance statistics
#' for two-way interactions and squares of covariates. Default is `FALSE`.
#' @param addlvariables additional variable for which balance statistics are to
#' be computed along with the covariates in the `matchit` object. Can be
#' entered in one of three ways: as a data frame of covariates with as many
#' rows as there were units in the original `matchit()` call, as a string
#' containing the names of variables in `data`, or as a right-sided
#' `formula` with the additional variables (and possibly their
#' transformations) found in `data`, the environment, or the
#' `matchit` object. Balance on squares and interactions of the additional
#' variables will be included if `interactions = TRUE`.
#' @param standardize `logical`; whether to compute standardized
#' (`TRUE`) or unstandardized (`FALSE`) statistics. The standardized
#' statistics are the standardized mean difference and the mean and maximum of
#' the difference in the (weighted) empirical cumulative distribution functions
#' (ECDFs). The unstandardized statistics are the raw mean difference and the
#' mean and maximum of the quantile-quantile (QQ) difference. Variance ratios
#' are produced either way. See Details below. Default is `TRUE`.
#' @param data a optional data frame containing variables named in
#' `addlvariables` if specified as a string or formula.
#' @param pair.dist `logical`; whether to compute average absolute pair
#' distances. For matching methods that don't include a `match.matrix`
#' component in the output (i.e., exact matching, coarsened exact matching,
#' full matching, and subclassification), computing pair differences can take a
#' long time, especially for large datasets and with many covariates. For other
#' methods (i.e., nearest neighbor, optimal, and genetic matching), computation
#' is fairly quick. Default is `FALSE` for subclassification and
#' `TRUE` otherwise.
#' @param un `logical`; whether to compute balance statistics for the
#' unmatched sample. Default `TRUE`; set to `FALSE` for more concise
#' output.
#' @param improvement `logical`; whether to compute the percent reduction
#' in imbalance. Default `FALSE`. Ignored if `un = FALSE`.
#' @param subclass after subclassification, whether to display balance for
#' individual subclasses, and, if so, for which ones. Can be `TRUE`
#' (display balance for all subclasses), `FALSE` (display balance only in
#' aggregate), or the indices (e.g., `1:6`) of the specific subclasses for
#' which to display balance. When anything other than `FALSE`, aggregate
#' balance statistics will not be displayed. Default is `FALSE`.
#' @param digits the number of digits to round balance statistics to.
#' @param x a `summay.matchit` or `summary.matchit.subclass` object;
#' the output of a call to `summary()`.
#' @param \dots ignored.
#'
#' @return For `matchit` objects, a `summary.matchit` object, which
#' is a list with the following components:
#'
#' \item{call}{the original call to [matchit()]}
#' \item{nn}{a matrix of the
#' sample sizes in the original (unmatched) and matched samples}
#' \item{sum.all}{if `un = TRUE`, a matrix of balance statistics for each
#' covariate in the original (unmatched) sample}
#' \item{sum.matched}{a matrix of
#' balance statistics for each covariate in the matched sample}
#' \item{reduction}{if `improvement = TRUE`, a matrix of the percent
#' reduction in imbalance for each covariate in the matched sample}
#'
#' For `match.subclass` objects, a `summary.matchit.subclass` object,
#' which is a list as above containing the following components:
#'
#' \item{call}{the original call to [matchit()]}
#' \item{sum.all}{if `un = TRUE`, a matrix of balance statistics for each covariate in the original
#' sample}
#' \item{sum.subclass}{if `subclass` is not `FALSE`, a list
#' of matrices of balance statistics for each subclass}
#' \item{sum.across}{a
#' matrix of balance statistics for each covariate computed using the
#' subclassification weights}
#' \item{reduction}{if `improvement = TRUE`, a
#' matrix of the percent reduction in imbalance for each covariate in the
#' matched sample}
#' \item{qn}{a matrix of sample sizes within each subclass}
#' \item{nn}{a matrix of the sample sizes in the original (unmatched) and
#' matched samples}
#'
#' @details
#' `summary()` computes a balance summary of a `matchit` object. This
#' include balance before and after matching or subclassification, as well as
#' the percent improvement in balance. The variables for which balance
#' statistics are computed are those included in the `formula`,
#' `exact`, and `mahvars` arguments to [matchit()], as well as the
#' distance measure if `distance` is was supplied as a numeric vector or
#' method of estimating propensity scores. The `X` component of the
#' `matchit` object is used to supply the covariates.
#'
#' The standardized mean differences are computed both before and after
#' matching or subclassification as the difference in treatment group means
#' divided by a standardization factor computed in the unmatched (original)
#' sample. The standardization factor depends on the argument supplied to
#' `estimand` in `matchit()`: for `"ATT"`, it is the standard
#' deviation in the treated group; for `"ATC"`, it is the standard
#' deviation in the control group; for `"ATE"`, it is the square root of
#' the average of the variances within each treatment group. The post-matching
#' mean difference is computed with weighted means in the treatment groups
#' using the matching or subclassification weights.
#'
#' The variance ratio is computed as the ratio of the treatment group
#' variances. Variance ratios are not computed for binary variables because
#' their variance is a function solely of their mean. After matching, weighted
#' variances are computed using the formula used in [cov.wt()]. The percent
#' reduction in bias is computed using the log of the variance ratios.
#'
#' The eCDF difference statistics are computed by creating a (weighted) eCDF
#' for each group and taking the difference between them for each covariate
#' value. The eCDF is a function that outputs the (weighted) proportion of
#' units with covariate values at or lower than the input value. The maximum
#' eCDF difference is the same thing as the Kolmogorov-Smirnov statistic. The
#' values are bounded at zero and one, with values closer to zero indicating
#' good overlap between the covariate distributions in the treated and control
#' groups. For binary variables, all eCDF differences are equal to the
#' (weighted) difference in proportion and are computed that way.
#'
#' The QQ difference statistics are computed by creating two samples of the
#' same size by interpolating the values of the larger one. The values are
#' arranged in order for each sample. The QQ difference for each quantile is
#' the difference between the observed covariate values at that quantile
#' between the two groups. The difference is on the scale of the original
#' covariate. Values close to zero indicate good overlap between the covariate
#' distributions in the treated and control groups. A weighted interpolation is
#' used for post-matching QQ differences. For binary variables, all QQ
#' differences are equal to the (weighted) difference in proportion and are
#' computed that way.
#'
#' The pair distance is the average of the absolute differences of a variable
#' between pairs. For example, if a treated unit was paired with four control
#' units, that set of units would contribute four absolute differences to the
#' average. Within a subclass, each combination of treated and control unit
#' forms a pair that contributes once to the average. The pair distance is
#' described in Stuart and Green (2008) and is the value that is minimized when
#' using optimal (full) matching. When `standardize = TRUE`, the
#' standardized versions of the variables are used, where the standardization
#' factor is as described above for the standardized mean differences. Pair
#' distances are not computed in the unmatched sample (because there are no
#' pairs). Because pair distance can take a while to compute, especially with
#' large datasets or for many covariates, setting `pair.dist = FALSE` is
#' one way to speed up `summary()`.
#'
#' The effective sample size (ESS) is a measure of the size of a hypothetical
#' unweighted sample with roughly the same precision as a weighted sample. When
#' non-uniform matching weights are computed (e.g., as a result of full
#' matching, matching with replacement, or subclassification), the ESS can be
#' used to quantify the potential precision remaining in the matched sample.
#' The ESS will always be less than or equal to the matched sample size,
#' reflecting the loss in precision due to using the weights. With non-uniform
#' weights, it is printed in the sample size table; otherwise, it is removed
#' because it does not contain additional information above the matched sample
#' size.
#'
#' After subclassification, the aggregate balance statistics are computed using
#' the subclassification weights rather than averaging across subclasses.
#'
#' All balance statistics (except pair differences) are computed incorporating
#' the sampling weights supplied to `matchit()`, if any. The unadjusted
#' balance statistics include the sampling weights and the adjusted balance
#' statistics use the matching weights multiplied by the sampling weights.
#'
#' When printing, `NA` values are replaced with periods (`.`), and
#' the pair distance column in the unmatched and percent balance improvement
#' components of the output are omitted.
#'
#' @seealso [summary()] for the generic method; [plot.summary.matchit()] for
#' making a Love plot from `summary()` output.
#'
#' \pkgfun{cobalt}{bal.tab.matchit}, which also displays balance for `matchit`
#' objects.
#'
#' @examples
#'
#' data("lalonde")
#' m.out <- matchit(treat ~ age + educ + married +
#'                    race + re74,
#'                  data = lalonde,
#'                  method = "nearest",
#'                  exact = ~ married,
#'                  replace = TRUE)
#'
#' summary(m.out, interactions = TRUE)
#'
#' s.out <- matchit(treat ~ age + educ + married +
#'                    race + nodegree + re74 + re75,
#'                  data = lalonde,
#'                  method = "subclass")
#'
#' summary(s.out, addlvariables = ~log(age) + I(re74==0))
#'
#' summary(s.out, subclass = TRUE)
#'

#' @exportS3Method summary matchit
summary.matchit <- function(object,
                            interactions = FALSE,
                            addlvariables = NULL,
                            standardize = TRUE,
                            data = NULL,
                            pair.dist = TRUE,
                            un = TRUE,
                            improvement = FALSE,
                            ...) {

  #Create covariate matrix; include caliper, exact, and mahvars
  X <- .process_X(object, addlvariables, data)

  treat <- object$treat
  weights <- object$weights
  s.weights <- {
    if (is_null(object$s.weights)) rep_with(1, weights)
    else object$s.weights
  }

  no_x <- is_null(X)

  if (no_x) {
    X <- matrix(1, nrow = length(treat), ncol = 1L,
                dimnames = list(names(treat), ".1"))
    nam <- colnames(X)
  }
  else {
    nam <- colnames(X)

    #Remove tics
    has_tics <- which(startsWith(nam, "`") & endsWith(nam, "`"))
    nam[has_tics] <- substr(nam[has_tics], 2, nchar(nam[has_tics]) - 1)
  }

  kk <- ncol(X)

  matched <- is_not_null(object$info$method)
  un <- un || !matched

  chk::chk_flag(interactions)
  chk::chk_flag(standardize)
  chk::chk_flag(pair.dist)
  chk::chk_flag(un)
  chk::chk_flag(improvement)

  s.d.denom <- {
    if (standardize) switch(object$estimand,
                            "ATT" = "treated",
                            "ATC" = "control",
                            "ATE" = "pooled")
    else NULL
  }

  ## Summary Stats
  if (un) {
    aa.all <- lapply(seq_len(kk), function(i) bal1var(X[,i], tt = treat, ww = NULL, s.weights = s.weights,
                                                      standardize = standardize, s.d.denom = s.d.denom))
    sum.all <- do.call("rbind", aa.all)
    dimnames(sum.all) <- list(nam, names(aa.all[[1L]]))

    if (no_x) sum.all <- sum.all[-1, , drop = FALSE]
    sum.all.int <- NULL
  }

  if (matched) {
    aa.matched <- lapply(seq_len(kk), function(i) bal1var(X[,i], tt = treat, ww = weights,
                                                          s.weights = s.weights,
                                                          subclass = object$subclass,
                                                          mm = object$match.matrix,
                                                          standardize = standardize,
                                                          s.d.denom = s.d.denom,
                                                          compute.pair.dist = pair.dist))
    sum.matched <- do.call("rbind", aa.matched)
    dimnames(sum.matched) <- list(nam, names(aa.matched[[1L]]))

    if (no_x) sum.matched <- sum.matched[-1, , drop = FALSE]
    sum.matched.int <- NULL
  }

  if (!no_x && interactions) {
    n.int <- kk * (kk + 1) / 2
    if (un) sum.all.int <- matrix(NA_real_, nrow = n.int, ncol = length(aa.all[[1L]]),
                                  dimnames = list(NULL, names(aa.all[[1]])))

    if (matched) sum.matched.int <- matrix(NA_real_, nrow = n.int, ncol = length(aa.matched[[1L]]),
                                           dimnames = list(NULL, names(aa.matched[[1L]])))

    to.remove <- rep.int(FALSE, n.int)
    int.names <- character(n.int)
    k <- 1
    for (i in seq_len(kk)) {
      for (j in i:kk) {
        x2 <- X[,i] * X[,j]
        if (all(abs(x2) < sqrt(.Machine$double.eps)) ||
            all(abs(x2 - X[,i]) < sqrt(.Machine$double.eps))) { #prevent interactions within same factors
          to.remove[k] <- TRUE
        }
        else {
          if (un) {
            sum.all.int[k,] <- bal1var(x2, tt = treat, ww = NULL, s.weights = s.weights,
                                       standardize = standardize, s.d.denom = s.d.denom)
          }
          if (matched) {
            sum.matched.int[k,] <- bal1var(x2, tt = treat, ww = weights, s.weights = s.weights,
                                           subclass = object$subclass, mm = object$match.matrix,
                                           standardize = standardize, s.d.denom = s.d.denom,
                                           compute.pair.dist = pair.dist)
          }
          if (i == j) {
            #Add superscript 2
            int.names[k] <- paste0(nam[i], "\u00B2")
          }
          else {
            int.names[k] <- paste(nam[i], nam[j], sep = " * ")
          }
        }
        k <- k + 1
      }
    }

    if (un) {
      rownames(sum.all.int) <- int.names
      sum.all <- rbind(sum.all, sum.all.int[!to.remove,,drop = FALSE])
    }

    if (matched) {
      rownames(sum.matched.int) <- int.names
      sum.matched <- rbind(sum.matched, sum.matched.int[!to.remove,,drop = FALSE])
    }
  }

  if (is_not_null(object$distance)) {
    if (un) {
      ad.all <- bal1var(object$distance, tt = treat, ww = NULL, s.weights = s.weights,
                        standardize = standardize, s.d.denom = s.d.denom)
      if (!exists("sum.all", inherits = FALSE)) {
        sum.all <- matrix(ad.all, nrow = 1L, dimnames = list("distance", names(ad.all)))
      }
      else {
        sum.all <- rbind(ad.all, sum.all)
        rownames(sum.all)[1L] <- "distance"
      }
    }
    if (matched) {
      ad.matched <- bal1var(object$distance, tt = treat, ww = weights, s.weights = s.weights,
                            subclass = object$subclass, mm = object$match.matrix, standardize = standardize,
                            s.d.denom = s.d.denom, compute.pair.dist = pair.dist)
      if (!exists("sum.matched", inherits = FALSE)) {
        sum.matched <- matrix(ad.matched, nrow = 1L,
                              dimnames = list("distance", names(ad.matched)))
      }
      else {
        sum.matched <- rbind(ad.matched, sum.matched)
        rownames(sum.matched)[1L] <- "distance"
      }
    }
  }

  ## Imbalance Reduction
  if (matched && un && improvement) {
    reduction <- matrix(NA_real_, nrow = nrow(sum.all), ncol = ncol(sum.all) - 2,
                        dimnames = list(rownames(sum.all), colnames(sum.all)[-(1:2)]))
    stat.all <- abs(sum.all[,-(1:2), drop = FALSE])
    stat.matched <- abs(sum.matched[,-(1:2), drop = FALSE])

    #Everything but variance ratios
    reduction[,-2] <- 100*(stat.all[,-2]-stat.matched[,-2])/stat.all[,-2]

    #Just variance ratios; turn to log first
    vr.all <- abs(log(stat.all[,2]))
    vr.matched <- abs(log(stat.matched[,2]))
    reduction[,2] <- 100*(vr.all-vr.matched)/vr.all

    reduction[stat.all == 0 & stat.matched == 0] <- 0
    reduction[stat.all == 0 & stat.matched > 0] <- -Inf
  }
  else {
    reduction <- NULL
  }

  #Sample size
  nn <- nn(treat, weights, object$discarded, s.weights)

  ## output
  res <- list(call = object$call,
              nn = nn,
              sum.all = if (un) sum.all,
              sum.matched = if (matched) sum.matched,
              reduction = reduction)

  class(res) <- "summary.matchit"

  res
}

#' @exportS3Method summary matchit.subclass
#' @rdname summary.matchit
summary.matchit.subclass <- function(object,
                                     interactions = FALSE,
                                     addlvariables = NULL,
                                     standardize = TRUE,
                                     data = NULL,
                                     pair.dist = FALSE,
                                     subclass = FALSE,
                                     un = TRUE,
                                     improvement = FALSE,
                                     ...) {

  #Create covariate matrix
  X <- .process_X(object, addlvariables, data)

  which.subclass <- subclass
  treat <- object$treat
  weights <- object$weights
  s.weights <- {
    if (is_null(object$s.weights)) rep_with(1, weights)
    else object$s.weights
  }
  subclass <- object$subclass

  nam <- colnames(X)

  kk <- ncol(X)
  subclasses <- levels(subclass)

  chk::chk_flag(interactions)
  chk::chk_flag(standardize)
  chk::chk_flag(pair.dist)
  chk::chk_flag(un)
  chk::chk_flag(improvement)

  s.d.denom <- {
    if (standardize) switch(object$estimand,
                            "ATT" = "treated",
                            "ATC" = "control",
                            "ATE" = "pooled")
    else NULL
  }

  if (isTRUE(which.subclass)) {
    which.subclass <- subclasses
  }
  else if (isFALSE(which.subclass)) {
    which.subclass <- NULL
  }
  else if (is.atomic(which.subclass) && all(which.subclass %in% seq_along(subclasses))) {
    which.subclass <- subclasses[which.subclass]
  }
  else {
    .err("`subclass` should be `TRUE`, `FALSE`, or a vector of subclass indices for which subclass balance is to be displayed")
  }

  matched <- TRUE #always compute aggregate balance so plot.summary can use it
  subs <- is_not_null(which.subclass)

  ## Aggregate Subclass
  #Use the estimated weights to compute aggregate balance.
  ## Summary Stats

  sum.all <- sum.matched <- sum.subclass <- reduction <- NULL

  if (un) {
    aa.all <- setNames(lapply(seq_len(kk), function(i) bal1var(X[,i], tt = treat, ww = NULL, s.weights = s.weights,
                                                               standardize = standardize, s.d.denom = s.d.denom)),
                       colnames(X))
    sum.all <- do.call("rbind", aa.all)
    dimnames(sum.all) <- list(nam, names(aa.all[[1L]]))

    sum.all.int <- NULL
  }

  if (matched) {
    aa.matched <- setNames(lapply(seq_len(kk), function(i) bal1var(X[,i], tt = treat, ww = weights, s.weights = s.weights,
                                                                   subclass = subclass, standardize = standardize,
                                                                   s.d.denom = s.d.denom, compute.pair.dist = pair.dist)),
                           colnames(X))
    sum.matched <- do.call("rbind", aa.matched)
    dimnames(sum.matched) <- list(nam, names(aa.matched[[1L]]))

    sum.matched.int <- NULL
  }

  if (interactions) {
    n.int <- kk * (kk + 1) / 2
    if (un) sum.all.int <- matrix(NA_real_, nrow = n.int, ncol = length(aa.all[[1L]]),
                                  dimnames = list(NULL, names(aa.all[[1L]])))
    if (matched) sum.matched.int <- matrix(NA_real_, nrow = n.int, ncol = length(aa.matched[[1L]]),
                                           dimnames = list(NULL, names(aa.matched[[1L]])))

    to.remove <- rep.int(FALSE, n.int)
    int.names <- character(n.int)
    k <- 1L
    for (i in seq_len(kk)) {
      for (j in i:kk) {
        x2 <- X[,i] * X[,j]
        if (all(abs(x2) < sqrt(.Machine$double.eps)) ||
            all(abs(x2 - X[,i]) < sqrt(.Machine$double.eps))) { #prevent interactions within same factors
          to.remove[k] <- TRUE
        }
        else {
          if (un) {
            sum.all.int[k,] <- bal1var(x2, tt = treat, ww = NULL, s.weights = s.weights,
                                       standardize = standardize, s.d.denom = s.d.denom)
          }
          if (matched) {
            sum.matched.int[k,] <- bal1var(x2, tt = treat, ww = weights, s.weights = s.weights,
                                           subclass = subclass, standardize = standardize,
                                           compute.pair.dist = pair.dist)
          }
          if (i == j) {
            int.names[k] <- paste0(nam[i], "\u00B2")
          }
          else {
            int.names[k] <- paste(nam[i], nam[j], sep = " * ")
          }
        }
        k <- k + 1L
      }
    }

    if (un) {
      rownames(sum.all.int) <- int.names
      sum.all <- rbind(sum.all, sum.all.int[!to.remove,,drop = FALSE])
    }
    if (matched) {
      rownames(sum.matched.int) <- int.names
      sum.matched <- rbind(sum.matched, sum.matched.int[!to.remove,,drop = FALSE])
    }
  }

  if (is_not_null(object$distance)) {
    if (un) {
      ad.all <- bal1var(object$distance, tt = treat, ww = NULL, s.weights = s.weights,
                        standardize = standardize, s.d.denom = s.d.denom)
      sum.all <- rbind(ad.all, sum.all)
      rownames(sum.all)[1L] <- "distance"
    }
    if (matched) {
      ad.matched <- bal1var(object$distance, tt = treat, ww = weights, s.weights = s.weights,
                            subclass = subclass, standardize = standardize,
                            s.d.denom = s.d.denom, compute.pair.dist = pair.dist)
      sum.matched <- rbind(ad.matched, sum.matched)
      rownames(sum.matched)[1L] <- "distance"
    }
  }

  ## Imbalance Reduction
  if (un && matched && improvement) {
    stat.all <- abs(sum.all[,-(1:2)])
    stat.matched <- abs(sum.matched[,-(1:2)])
    reduction <- 100 * (stat.all - stat.matched) / stat.all

    reduction[stat.all == 0 & stat.matched == 0] <- 0
    reduction[stat.all == 0 & stat.matched > 0] <- -Inf
  }

  ## By Subclass
  if (subs) {
    sum.subclass <- lapply(which.subclass, function(s) {

      #bal1var.subclass only returns unmatched stats, which is all we need within
      #subclasses. Otherwise, identical to matched stats.
      aa <- setNames(lapply(seq_len(kk), function(i) {
        bal1var.subclass(X[,i], tt = treat, s.weights = s.weights,
                         subclass = subclass, s.d.denom = s.d.denom,
                         standardize = standardize, which.subclass = s)
      }), colnames(X))

      sum.sub <- matrix(NA_real_, nrow = kk, ncol = ncol(aa[[1L]]), dimnames = list(nam, colnames(aa[[1L]])))

      sum.sub.int <- NULL
      for (i in seq_len(kk)) {
        sum.sub[i,] <- aa[[i]]
      }
      if (interactions) {
        sum.sub.int <- matrix(NA_real_, nrow = kk * (kk + 1) / 2, ncol = length(aa[[1L]]),
                              dimnames = list(NULL, names(aa[[1L]])))
        to.remove <- rep.int(FALSE, nrow(sum.sub.int))
        int.names <- character(nrow(sum.sub.int))
        k <- 1L
        for (i in seq_len(kk)) {
          for (j in i:kk) {
            if (!to.remove[k]) { #to.remove defined above
              x2 <- X[,i] * X[,j]
              jqoi <- bal1var.subclass(x2, tt = treat, s.weights = s.weights,
                                       subclass = subclass, s.d.denom = s.d.denom,
                                       standardize = standardize, which.subclass = s)
              sum.sub.int[k,] <- jqoi
              if (i == j) {
                int.names[k] <- paste0(nam[i], "\u00B2")
              }
              else {
                int.names[k] <- paste(nam[i], nam[j], sep = " * ")
              }
            }
            k <- k + 1
          }
        }
        rownames(sum.sub.int) <- int.names

        sum.sub <- rbind(sum.sub, sum.sub.int[!to.remove,,drop = FALSE])
      }

      if (is_not_null(object$distance)) {
        ad <- bal1var.subclass(object$distance, tt = treat, s.weights = s.weights, subclass = subclass,
                               s.d.denom = s.d.denom, standardize = standardize, which.subclass = s)
        sum.sub <- rbind(ad, sum.sub)
        rownames(sum.sub)[1L] <- "distance"
      }

      sum.sub
    })
    names(sum.subclass) <- paste("Subclass", which.subclass)
  }

  ## Sample size
  qn <- qn(treat, subclass, object$discarded)
  nn <- nn(treat, weights, object$discarded, s.weights)

  ## output
  res <- list(call = object$call,
              sum.all = sum.all,
              sum.across = sum.matched,
              sum.subclass = sum.subclass,
              reduction = reduction,
              qn = qn,
              nn = nn)

  class(res) <- c("summary.matchit.subclass", "summary.matchit")

  res
}

#' @exportS3Method print summary.matchit
#' @rdname summary.matchit
print.summary.matchit <- function(x, digits = max(3, getOption("digits") - 3),
                                  ...) {

  if (is_not_null(x$call)) {
    cat("\nCall:", deparse(x$call), sep = "\n")
  }

  if (is_not_null(x$sum.all)) {
    cat("\nSummary of Balance for All Data:\n")
    print(round_df_char(x$sum.all[,-7, drop = FALSE], digits, pad = "0", na_vals = "."),
          right = TRUE, quote = FALSE)
  }

  if (is_not_null(x$sum.matched)) {
    cat("\nSummary of Balance for Matched Data:\n")
    if (all(is.na(x$sum.matched[,7]))) x$sum.matched <- x$sum.matched[,-7,drop = FALSE] #Remove pair dist if empty
    print(round_df_char(x$sum.matched, digits, pad = "0", na_vals = "."),
          right = TRUE, quote = FALSE)
  }
  if (is_not_null(x$reduction)) {
    cat("\nPercent Balance Improvement:\n")
    print(round_df_char(x$reduction[,-5, drop = FALSE], 1, pad = "0", na_vals = "."), right = TRUE,
          quote = FALSE)
  }
  if (is_not_null(x$nn)) {
    cat("\nSample Sizes:\n")
    nn <- x$nn
    if (isTRUE(all.equal(nn["All (ESS)",], nn["All",]))) {
      #Don't print ESS if same as full SS
      nn <- nn[rownames(nn) != "All (ESS)",,drop = FALSE]
    }
    if (isTRUE(all.equal(nn["Matched (ESS)",], nn["Matched",]))) {
      #Don't print ESS if same as matched SS
      nn <- nn[rownames(nn) != "Matched (ESS)",,drop = FALSE]
    }
    print(round_df_char(nn, 2, pad = " ", na_vals = "."), right = TRUE,
          quote = FALSE)
  }
  cat("\n")
  invisible(x)
}

#' @exportS3Method print summary.matchit.subclass
print.summary.matchit.subclass <- function(x, digits = max(3, getOption("digits") -  3), ...){

  if (is_not_null(x$call)) {
    cat("\nCall:", deparse(x$call), sep = "\n")
  }

  if (is_not_null(x$sum.all)) {
    cat("\nSummary of Balance for All Data:\n")
    print(round_df_char(x$sum.all[,-7, drop = FALSE], digits, pad = "0", na_vals = "."),
          right = TRUE, quote = FALSE)
  }

  if (is_not_null(x$sum.subclass)) {
    cat("\nSummary of Balance by Subclass:\n")
    for (s in seq_along(x$sum.subclass)) {
      cat(paste0("\n- ", names(x$sum.subclass)[s], "\n"))
      print(round_df_char(x$sum.subclass[[s]][,-7, drop = FALSE], digits, pad = "0", na_vals = "."),
            right = TRUE, quote = FALSE)
    }
    if (is_not_null(x$qn)) {
      cat("\nSample Sizes by Subclass:\n")
      print(round_df_char(x$qn, 2, pad = " ", na_vals = "."),
            right = TRUE, quote = FALSE)
    }
  }
  else {
    if (is_not_null(x$sum.across)) {
      cat("\nSummary of Balance Across Subclasses\n")
      if (all(is.na(x$sum.across[,7]))) x$sum.across <- x$sum.across[,-7,drop = FALSE]
      print(round_df_char(x$sum.across, digits, pad = "0", na_vals = "."),
            right = TRUE, quote = FALSE)
    }
    if (is_not_null(x$reduction)) {
      cat("\nPercent Balance Improvement:\n")
      print(round_df_char(x$reduction[,-5, drop = FALSE], 1, pad = "0", na_vals = "."),
            right = TRUE, quote = FALSE)
    }

    if (is_not_null(x$nn)) {
      cat("\nSample Sizes:\n")
      nn <- x$nn
      if (isTRUE(all.equal(nn["All (ESS)",], nn["All",]))) {
        #Don't print ESS if same as full SS
        nn <- nn[rownames(nn) != "All (ESS)",,drop = FALSE]
      }
      if (isTRUE(all.equal(nn["Matched (ESS)",], nn["Matched",]))) {
        #Don't print ESS if same as matched SS
        nn <- nn[rownames(nn) != "Matched (ESS)",,drop = FALSE]
      }
      print(round_df_char(nn, 2, pad = " ", na_vals = "."),
            right = TRUE, quote = FALSE)
    }
  }
  cat("\n")
}

.process_X <- function(object, addlvariables = NULL, data = NULL) {

  X <- {
    if (is_null(object$X)) matrix(nrow = length(object$treat), ncol = 0)
    else get_covs_matrix(data = object$X)
  }

  if (is_null(addlvariables)) {
    return(X)
  }

  #Attempt to extract data from matchit object; same as match_data()
  data.found <- FALSE
  for (i in 1:4) {
    if (i == 2L) {
      data <- try(eval(object$call$data, envir = environment(object$formula)), silent = TRUE)
    }
    else if (i == 3L) {
      data <- try(eval(object$call$data, envir = parent.frame()), silent = TRUE)
    }
    else if (i == 4L) {
      data <- object[["model"]][["data"]]
    }

    if (!null_or_error(data) && length(dim(data)) == 2L && nrow(data) == length(object[["treat"]])) {
      data.found <- TRUE
      break
    }
  }

  if (is.character(addlvariables)) {
    if (is_null(data) || !is.data.frame(data)) {
      .err("if `addlvariables` is specified as a string, a data frame argument must be supplied to `data`")
    }

    if (!all(hasName(data, addlvariables))) {
      .err("all variables in `addlvariables` must be in `data`")
    }

    addlvariables <- data[addlvariables]
  }
  else if (rlang::is_formula(addlvariables)) {
    if (is_not_null(data) && is.data.frame(data)) {
      vars.in.formula <- all.vars(addlvariables)
      data <- cbind(data[names(data) %in% vars.in.formula],
                    object$X[names(object$X) %in% setdiff(vars.in.formula, names(data))])
    }
    else {
      data <- object$X
    }
  }
  else if (!is.matrix(addlvariables) && !is.data.frame(addlvariables)) {
    .err("the argument to `addlvariables` must be in one of the accepted forms. See `?summary.matchit` for details")
  }


  if (af <- rlang::is_formula(addlvariables)) {
    addvariables_f <- addlvariables
    addlvariables <- model.frame(addvariables_f, data = data, na.action = "na.pass")
  }

  if (nrow(addlvariables) != length(object$treat)) {
    if (is_null(data) || data.found) {
      .err("variables specified in `addlvariables` must have the same number of units as are present in the original call to `matchit()`")
    }
    else {
      .err("`data` must have the same number of units as are present in the original call to `matchit()`")
    }
  }

  k <- ncol(addlvariables)
  for (i in seq_len(k)) {
    if (anyNA(addlvariables[[i]]) || (is.numeric(addlvariables[[i]]) && any(!is.finite(addlvariables[[i]])))) {
      covariates.with.missingness <- names(addlvariables)[i:k][vapply(i:k, function(j) anyNA(addlvariables[[j]]) ||
                                                                        (is.numeric(addlvariables[[j]]) &&
                                                                           any(!is.finite(addlvariables[[j]]))),
                                                                      logical(1L))]
      .err(paste0("Missing and non-finite values are not allowed in `addlvariables`. Variables with missingness or non-finite values:\n\t",
                  paste(covariates.with.missingness, collapse = ", ")), tidy = FALSE)
    }

    if (is.character(addlvariables[[i]])) {
      addlvariables[[i]] <- factor(addlvariables[[i]])
    }
  }

  if (af) {
    addlvariables <- get_covs_matrix(addvariables_f, data = data)
  }
  else {
    addlvariables <- get_covs_matrix(data = addlvariables)
  }

  # addl_assign <- get_assign(addlvariables)
  cbind(X, addlvariables[, setdiff(colnames(addlvariables), colnames(X)), drop = FALSE])

}
