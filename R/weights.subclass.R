weights.subclass <-function(psclass, treat) {

  n <- length(treat)
  labels <- names(treat)
  tlabels <- labels[treat==1]
  clabels <- labels[treat==0]

  weights <- rep(0,length(treat))
  names(weights) <- labels
  weights[tlabels][!is.na(psclass[tlabels])] <- 1

  classes <- unique(psclass)
  for(j in classes){
        qn0 <- sum(treat==0 & psclass==j)
        qn1 <- sum(treat==1 & psclass==j)
        weights[treat==0 & psclass==j] <- qn1/qn0
      }
    if (sum(weights[treat==0])==0) {
      weights[treat==0] <- rep(0, length(weights[clabels]))
    } else {
      # Number of C units that were matched to at least 1 T
      num.cs <- sum(weights[clabels]>0)
      weights[clabels] <- weights[clabels]*num.cs/sum(weights[clabels])
    }

  if (sum(weights)==0) {
    stop("Error: No units were matched")
  } else if (sum(weights[tlabels])==0){
    stop("No treated units were matched",call.=FALSE)
  } else if (sum(weights[clabels])==0){
    stop("No control units were matched",call.=FALSE)
  }
  cat("Done\n")
  list(psweights=weights)
}

