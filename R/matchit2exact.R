matchit2exact <- function(treat, X, data, pscore, discarded, verbose=FALSE, ...){

   if(verbose)
    cat("Exact matching... \n")

  n <- length(treat)
  xx <- apply(X, 1, function(x) paste(x, collapse = "\r"))
  xx1 <- xx[treat==1]
  xx0 <- xx[treat==0]
  cc <- unique(xx1)
  cc <- cc[cc%in%xx0]
  ncc <- length(cc)
  
  psclass <- rep(NA,n)
  names(psclass) <- names(treat)
  for(i in 1:ncc)
    psclass[xx==cc[i]] <- i

  res <- list(subclass = psclass, weights = weights.subclass(psclass, treat))
  class(res) <- c("matchit.exact", "matchit")
  return(res)
}
