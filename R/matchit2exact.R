matchit2exact <- function(treat, X, dis,  ...){

  n <- length(treat)
  
  xx <- apply(X, 1, function(x) paste(x, collapse = "\r"))

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

  names(psclass) <- names(treat)
  
  psweights <- weights.subclass(psclass=psclass, treat=treat)$psweights
 
  res <- list(psclass=psclass,  q.cut=NULL, psweights=psweights)
  class(res) <- c("matchit.subclass", "matchit")
  return(res)
}
