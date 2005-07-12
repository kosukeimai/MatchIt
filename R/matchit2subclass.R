matchit2subclass <- function(treat, X, data, pscore, discarded,
                             match.matrix=NULL, subclass=6, sub.by="treat",
                             verbose = FALSE){
  in.sample <- !discarded
  if(verbose)
    cat("Subclassifying...")  
  n <- length(treat)
  ## Matching & Subclassification
  if(!is.null(match.matrix)){
    #match.matrix <- match.matrix[in.sample[treat==1],,drop=F]
    t.units <- row.names(match.matrix)[in.sample[treat==1]==1]
    c.units <- na.omit(as.vector(as.matrix(match.matrix)))
    matched <-c(t.units,c.units)
    matched <- names(treat)%in%matched
  } else
    matched <- rep(TRUE,n)
  names(matched) <- names(treat)
  m1 <- matched[treat==1]
  m0 <- matched[treat==0]
  p1 <- pscore[treat==1][m1]
  p0 <- pscore[treat==0][m0]
  ## Settting Cut Points
  if(length(subclass)!=1 | (length(subclass)==1 & all(subclass<1))) {
    subclass <- sort(subclass)
    if (subclass[1]==0)
      subclass <- subclass[-1]
    if (subclass[length(subclass)]==1)
      subclass <- subclass[-length(subclass)]
    if(sub.by=="treat") 
      q <- c(0,quantile(p1,probs=c(subclass)),1)
    else if(sub.by=="control") 
      q <- c(0,quantile(p0,probs=c(subclass)),1)
    else if(sub.by=="all") 
      q <- c(0,quantile(pscore,probs=c(subclass)),1)
    else 
      stop("Invalid input for sub.by")
  }
  else {
    if(subclass<=0){stop("Subclass must be a positive vector",call.=FALSE)}
    sprobs <- seq(0,1,length=(round(subclass)+1))
    sprobs <- sprobs[2:(length(sprobs)-1)]
    if(sub.by=="treat")   
      q <- c(0,quantile(p1,probs=sprobs,na.rm=TRUE),1)
    else if(sub.by=="control") 
      q <- c(0,quantile(p0,probs=sprobs,na.rm=TRUE),1)
    else if(sub.by=="all") 
      q <- c(0,quantile(pscore,probs=sprobs,na.rm=TRUE),1)
    else 
      stop("Must specify a valid sub.by",call.=FALSE)
  }
  ## Calculating Subclasses
  qbins <- length(q)-1
  psclass <- rep(0,n)
  names(psclass) <- names(treat)
  for (i in 1:qbins){
    q1 <- q[i]
    q2 <- q[i+1]
    psclass <- psclass+i*as.numeric(pscore<q2 & pscore>=q1)
  }
  
  ## No subclass for discarded or unmatched units
  psclass[in.sample==0] <- NA
  psclass[!matched] <- NA
  if(verbose){cat("Done\n")}
  res <- list(subclass = psclass, q.cut = q,
              weights = weights.subclass(psclass, treat))

  class(res) <- c("matchit.subclass", "matchit")
  return(res)
}
