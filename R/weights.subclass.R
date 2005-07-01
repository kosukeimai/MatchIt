weights.subclass <-function(psclass, treat) {

  print(psclass)

  n <- length(treat)
  labels <- names(treat)
  tlabels <- labels[treat==1]
  clabels <- labels[treat==0]

  weights <- rep(0,length(treat))
  names(weights) <- labels
  weights[tlabels][!is.na(psclass[tlabels])] <- 1
 
  classes <- unique(psclass)
  for(j in classes){
        qn0 <- sum(treat==0 & psclass==j, na.rm=T)
        qn1 <- sum(treat==1 & psclass==j, na.rm=T)
        weights[treat==0 & psclass==j] <- qn1/qn0
        print(c(j, qn0, qn1))
      }
    if (sum(weights[treat==0], na.rm=T)==0) {
      weights[treat==0] <- rep(0, length(weights[clabels]))
    } else {
      # Number of C units that were matched to at least 1 T
      num.cs <- sum(weights[clabels]>0, na.rm=T)
      weights[clabels][!is.na(psclass[clabels])] <- 
	weights[clabels][!is.na(psclass[clabels])]*num.cs/sum(weights[clabels], na.rm=T)
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

