discard <- function(treat, pscore, option, X) {

  n.obs <- length(treat)
  pmax0 <- max(pscore[treat==0])
  pmax1 <- max(pscore[treat==1])
  pmin0 <- min(pscore[treat==0])
  pmin1 <- min(pscore[treat==1])
  if (option == "none")         # keep all units
    discarded <- rep(FALSE, n.obs)
  else if (option == "both")    # discard units outside of common support
    discarded <- (pscore < max(pmin0, pmin1) | pscore > min(pmax0, pmax1))
  else if (option == "control") # discard control units only
    discarded <- (pscore < pmin1 | pscore > pmax1)
  else if (option == "treat")   # discard treated units only
    discarded <- (pscore < pmin0 | pscore > pmax0)
  else if (option == "convex.hull"){ # discard units not in T convex hull
    if (!("whatif" %in% .packages(all = TRUE)))
      install.packages("whatif", CRAN="http://gking.harvard.edu")
    if (!("lpSolve" %in% .packages(all = TRUE)))
      install.packages("lpSolve")
    require(whatif)
    require(lpSolve)
    wif <- whatif(cfact = X[treat==0,], data = X[treat==1,])
    print(wif$in.hull)
    discarded <- rep(FALSE, n.obs)
    discarded[treat==0] <- !wif$in.hull
  }  else 
    stop("invalid input for `discard'")
  names(discarded) <- names(treat)
  return(discarded)
}
