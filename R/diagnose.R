diagnose <-function(formula, match.matrix, pscore, in.sample, data,
                    exact = FALSE, mahvars = NULL, subclass = 0,
                    psclass = NULL, nearest = TRUE, q.cut = NULL,
                    counter = TRUE){
  if(counter){cat("Calculating summary statistics...")}
  treata <- model.frame(formula,data)[,1,drop=FALSE]
  treat <- as.vector(treata[,1])
  names(treat) <- row.names(treata)
  n <- length(treat)
  n0 <- length(treat[treat==0])
  n1 <- length(treat[treat==1])

  if(is.null(names(treat))){names(treat) <- seq(1,n)}
  labels <- names(treat)
  clabels <- labels[treat==0]
  tlabels <- labels[treat==1]

  # setting up weights
  if(nearest==TRUE){
    match.matrix <- match.matrix[in.sample[treat==1],,drop=F]
    num.matches <- dim(match.matrix)[2]-apply(as.matrix(match.matrix), 1, function(x) { sum(is.na(x)) })
    names(num.matches) <- tlabels[in.sample[treat==1]]

    t.units <- row.names(match.matrix)[num.matches>0]
    c.units <- na.omit(as.vector(as.matrix(match.matrix)))
    matched <-c(t.units,unique(c.units))
    
    weights <- rep(0,length(treat))
    names(weights) <- labels
    weights[t.units] <- 1
    
    for (cont in clabels) {
      treats <- na.omit(row.names(match.matrix)[cont==match.matrix[,1]])
      if (dim(match.matrix)[2]>1) {
        for (j in 2:dim(match.matrix)[2]) 
          treats <- c(na.omit(row.names(match.matrix)[cont==match.matrix[,j]]),treats)	
      }
      weights[cont] <- length(unique(treats))
    }
    
    if (sum(weights[clabels])==0){
      weights[clabels] <- rep(0, length(weights[clabels]))
    } else {
      weights[clabels] <-length(weights[clabels][weights[clabels]>0])*weights[clabels]/sum(weights[clabels][weights[clabels]>0])
    }
    
  } else if(identical(exact,TRUE)) {
    
    weights <- rep(0,length(treat))
    names(weights) <- labels
    t.units <- names(psclass)[psclass>0 & treat==1]
    weights[t.units] <- 1
    for(j in 1:max(psclass)){
        qn0 <- sum(treat==0 & psclass==j)
        qn1 <- sum(treat==1 & psclass==j)
        weights[treat==0 & psclass==j] <- qn1/qn0
      }
    if (sum(weights[treat==0])==0) {
      weights[treat==0] <- rep(0, length(weights[clabels]))
    } else {
      weights[clabels] <- length(weights[treat==0][weights[treat==0]>0]) *
        weights[treat==0]/sum(weights[treat==0][weights[treat==0]>0])
    }
    matched <- names(treat)
    
  } else {
    matched <- names(treat)
    weights <- rep(1,length(treat))
    names(weights) <- names(treat)
  }
  weights[!in.sample] <- 0
  if (sum(weights)==0) {
    stop("Error: No units were matched")
  } else if (sum(weights[tlabels])==0){
    stop("No treated units were matched",call.=FALSE)
  } else if (sum(weights[clabels])==0){
    stop("No control units were matched",call.=FALSE)
  }
  cat("Done\n")
  list(weights=weights)
}

