matchit2exact <- function(formula, data, counter=TRUE, ...){

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
  
   z <- list(match.matrix = NULL, in.sample = NULL, pscore=NULL, assign.model=NULL, psclass=psclass)
  class(z) <- "matchit"
  return(z)
}
