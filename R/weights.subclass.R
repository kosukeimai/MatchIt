weights.subclass <- function(psclass, treat) {

  ttt <- treat[!is.na(psclass)]
  classes <- na.omit(psclass)

  # Determine which number code corresponds to the treated group
  if(is.factor(ttt)){
      treated <- as.numeric(ttt) == 2
  } else {
      treated <- ttt == 1
  }
    
  n <- length(ttt)
  labels <- names(ttt)
  tlabels <- labels[treated]
  clabels <- labels[!treated]
  
  weights <- rep(0, n)
  names(weights) <- labels
  weights[tlabels] <- 1
 
  for(j in unique(classes)){
    qn0 <- sum(!treated & classes==j)
    qn1 <- sum(treated & classes==j)
    weights[!treated & classes==j] <- qn1/qn0
  }
  if (sum(weights[!treated])==0)
    weights[!treated] <- rep(0, length(weights[clabels]))
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

