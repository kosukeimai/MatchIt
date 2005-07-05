weights.subclass <-function(psclass, treat) {

  n <- length(treat)
  labels <- names(treat)
  tlabels <- labels[treat==1]
  clabels <- labels[treat==0]

  weights <- rep(0,length(treat))
  names(weights) <- labels
  weights[tlabels][psclass[tlabels]!=0] <- 1
 
  classes <- unique(psclass)
  for(j in classes){
    qn0 <- sum(treat==0 & psclass==j)
    qn1 <- sum(treat==1 & psclass==j)
    weights[treat==0 & psclass==j] <- qn1/qn0
  }
  if (sum(weights[treat==0], na.rm=T)==0)
    weights[treat==0] <- rep(0, length(weights[clabels]))
  else {
    ## Number of C units that were matched to at least 1 T
    num.cs <- sum(weights[clabels][psclass[clabels]!=0]>0)
    weights[clabels][psclass[clabels]!=0] <- 
      weights[clabels][psclass[clabels]!=0]*num.cs/sum(weights[clabels][psclass[clabels]!=0])
  }
  
  if (sum(weights)==0) 
    stop("No units were matched")
  else if (sum(weights[tlabels])==0)
    stop("No treated units were matched")
  else if (sum(weights[clabels])==0)
    stop("No control units were matched")
  
  return(list(psweights=weights))
}

