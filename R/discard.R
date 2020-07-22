discard <- function(treat, pscore = NULL, option = NULL, covs = NULL) {

  n.obs <- length(treat)

  if (is.logical(option) && length(option) == n.obs) {
    # user input
    return(option)
  }
  else if (length(option) > 1 || !is.character(option)) {
    stop("Invalid input for discard.", call. = FALSE)
  }
  else if (length(option) == 0 || option == "none"){
    # keep all units
    discarded <- rep(FALSE, n.obs)
  }
  else if (option %in% c("both", "control", "treat")) {
    pmax0 <- max(pscore[treat==0])
    pmax1 <- max(pscore[treat==1])
    pmin0 <- min(pscore[treat==0])
    pmin1 <- min(pscore[treat==1])

    if (option == "both")    # discard units outside of common support
      discarded <- (pscore < max(pmin0, pmin1) | pscore > min(pmax0, pmax1))
    else if (option == "control") # discard control units only
      discarded <- (pscore < pmin1 | pscore > pmax1)
    else if (option == "treat")   # discard treated units only
      discarded <- (pscore < pmin0 | pscore > pmax0)
  }
  else if (option %in% c("hull.control", "hull.treat", "hull.both")) {
    ## convex hull stuff
    check.package("WhatIf")
    X <- model.matrix(reformulate(names(covs), intercept = FALSE), data = covs,
                      contrasts.arg = lapply(Filter(is.factor, covs),
                                             function(x) contrasts(x, contrasts = nlevels(x) == 1)))
    discarded <- rep(FALSE, n.obs)
    if (option == "hull.control"){ # discard units not in T convex hull
      wif <- WhatIf::whatif(cfact = X[treat==0,], data = X[treat==1,])
      discarded[treat==0] <- !wif$in.hull
    } else if (option == "hull.treat") {
      wif <- WhatIf::whatif(cfact = X[treat==1,], data = X[treat==0,])
      discarded[treat==1] <- !wif$in.hull
    } else if (option == "hull.both"){ # discard units not in T&C convex hull
      wif <- WhatIf::whatif(cfact = cbind(1-treat, X), data = cbind(treat, X))
      discarded <- !wif$in.hull
    }
    else
      stop("invalid input for discard")
  } else
  stop("Invalid input for discard.", call. = FALSE)

  names(discarded) <- names(treat)

  return(discarded)
}
