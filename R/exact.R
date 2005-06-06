exact <- function(formula, data, counter=TRUE, ...){

  if(counter){cat("Exact matching...")}
  data <- eval(data,parent.frame())
  treata <- model.frame(formula,data)[,1,drop=FALSE]
  treat <- as.vector(treata[,1])
  names(treat) <- row.names(treata)
  covariates <- model.frame(delete.response(terms(formula)),data)[,,drop=FALSE]
  n <- length(treat)
  n1 <- length(treat[treat==1])
  n0 <- length(treat[treat==0])
  covariates <- as.data.frame(covariates)
  xx <- apply(covariates, 1, function(x) paste(x, collapse = "\r"))
  xx1 <- xx[treat==1]
  xx0 <- xx[treat==0]
  cc <- unique(xx1)
  cc <- cc[cc%in%xx0]
  ncc <- length(cc)
  psclass <- rep(0,n)
  names(psclass) <- names(treat)
  for(i in 1:ncc){
    psclass[xx==cc[i]] <- i
  }
  if(counter){cat("Done\n")}

  ##########################################################################
  ## Calculate weights
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

   z <- list(match.matrix = NULL, in.sample = NULL, pscore=NULL, assign.model=NULL, psclass=psclass, 
		psweights=weights)
  class(z) <- "matchit"
  return(z)
}
