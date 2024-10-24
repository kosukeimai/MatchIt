discard <- function(treat, pscore = NULL, option = NULL) {

  n.obs <- length(treat)

  if (is_null(option)){
    # keep all units
    return(setNames(rep(FALSE, n.obs), names(treat)))
  }

  if (is.logical(option) && length(option) == n.obs && !anyNA(option)) {
    # user input
    return(setNames(option, names(treat)))
  }

  if (!chk::vld_string(option)) {
    .err('`discard` must be "none", "both", "control", "treated" or a logical vector of observations to discard')
  }

  option <- match_arg(option, c("none", "both", "control", "treated"))

  if (option == "none") {
    # keep all units
    return(setNames(rep(FALSE, n.obs), names(treat)))
  }

  if (is_null(pscore)) {
    .err('`discard` must be a logical vector or "none" in the absence of a propensity score')
  }

  if (is.matrix(pscore)) {
    .err('`discard` must be a logical vector or "none" when `distance` is supplied as a matrix')
  }

  pmax0 <- max(pscore[treat==0])
  pmax1 <- max(pscore[treat==1])
  pmin0 <- min(pscore[treat==0])
  pmin1 <- min(pscore[treat==1])

  if (option == "both")    # discard units outside of common support
    discarded <- (pscore < max(pmin0, pmin1) | pscore > min(pmax0, pmax1))
  else if (option == "control") # discard control units only
    discarded <- (pscore < pmin1 | pscore > pmax1)
  else if (option == "treated")   # discard treated units only
    discarded <- (pscore < pmin0 | pscore > pmax0)

  # NOTE: WhatIf package has been removed from CRAN, so hull options won't work
  # else if (option %in% c("hull.control", "hull.treat", "hull.both")) {
  #   ## convex hull stuff
  #   check.package("WhatIf")
  #   X <- model.matrix(reformulate(names(covs), intercept = FALSE), data = covs,
  #                     contrasts.arg = lapply(Filter(is.factor, covs),
  #                                            function(x) contrasts(x, contrasts = nlevels(x) == 1)))
  #   discarded <- rep(FALSE, n.obs)
  #   if (option == "hull.control"){ # discard units not in T convex hull
  #     wif <- WhatIf::whatif(cfact = X[treat==0,], data = X[treat==1,])
  #     discarded[treat==0] <- !wif$in.hull
  #   } else if (option == "hull.treat") {
  #     wif <- WhatIf::whatif(cfact = X[treat==1,], data = X[treat==0,])
  #     discarded[treat==1] <- !wif$in.hull
  #   } else if (option == "hull.both"){ # discard units not in T&C convex hull
  #     wif <- WhatIf::whatif(cfact = cbind(1-treat, X), data = cbind(treat, X))
  #     discarded <- !wif$in.hull
  #   }
  # }

  setNames(discarded, names(treat))
}
