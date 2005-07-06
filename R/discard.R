discard <- function(treat, pscore, option) {

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
  else
    stop("invalid input for `discard'")
  
  return(discarded)
}
