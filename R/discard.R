discard <- function(treat, pscore, option, X) {

  n.obs <- length(treat)
  pmax0 <- max(pscore[treat==0])
  pmax1 <- max(pscore[treat==1])
  pmin0 <- min(pscore[treat==0])
  pmin1 <- min(pscore[treat==1])
  if (is.logical(option))       # user input
    return(option)
  else if (option == "none")    # keep all units
    discarded <- rep(FALSE, n.obs)
  else if (option == "both")    # discard units outside of common support
    discarded <- (pscore < max(pmin0, pmin1) | pscore > min(pmax0, pmax1))
  else if (option == "control") # discard control units only
    discarded <- (pscore < pmin1 | pscore > pmax1)
  else if (option == "treat")   # discard treated units only
    discarded <- (pscore < pmin0 | pscore > pmax0)
  else if (any(grep(option, c("hull.control", "hull.treat", "hull.both")))) {
    ## convext hull stuff
    if (!("WhatIf" %in% .packages(all = TRUE)))
      install.packages("WhatIf")
    if (!("lpSolve" %in% .packages(all = TRUE)))
      install.packages("lpSolve")
    require(WhatIf)
    require(lpSolve)
    discarded <- rep(FALSE, n.obs)
    if (option == "hull.control"){ # discard units not in T convex hull
      wif <- whatif(cfact = X[treat==0,], data = X[treat==1,])
      discarded[treat==0] <- !wif$in.hull
    } else if (option == "hull.treat") {
      wif <- whatif(cfact = X[treat==1,], data = X[treat==0,])
      discarded[treat==1] <- !wif$in.hull
    } else if (option == "hull.both"){ # discard units not in T&C convex hull
      wif <- whatif(cfact = cbind(1-treat, X), data = cbind(treat, X))
      discarded <- !wif$in.hull
    }
    else
      stop("invalid input for `discard'")
  } else 
  stop("invalid input for `discard'")
  names(discarded) <- names(treat)
  return(discarded)
}
