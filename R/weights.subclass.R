weights.subclass <- function(psclass, treat, estimand = "ATT") {

  ttt <- treat[!is.na(psclass)]
  classes <- psclass[!is.na(psclass)]

  n <- length(ttt)
  labels <- names(ttt)
  tlabels <- labels[ttt==1]
  clabels <- labels[ttt==0]

  weights <- rep(0, n)
  names(weights) <- labels
  weights[tlabels] <- 1

  for (j in unique(classes)){
    qn0 <- sum(ttt==0 & classes==j)
    qn1 <- sum(ttt==1 & classes==j)
    weights[ttt==0 & classes==j] <- qn1/qn0
  }
  if (sum(weights[ttt==0])==0)
    weights[ttt==0] <- rep(0, length(weights[clabels]))
  else {
    ## Number of C units that were matched to at least 1 T
    num.cs <- sum(weights[clabels] > 0)
    weights[clabels] <- weights[clabels]*num.cs/sum(weights[clabels])
  }

  if (any(is.na(psclass))) {
    tmp <- rep(0, sum(is.na(psclass)))
    names(tmp) <- names(treat[is.na(psclass)])
    weights <- c(weights, tmp)[names(treat)]
  }

  if (sum(weights)==0)
    stop("No units were matched")
  else if (sum(weights[tlabels])==0)
    stop("No treated units were matched")
  else if (sum(weights[clabels])==0)
    stop("No control units were matched")

  return(weights)
}

