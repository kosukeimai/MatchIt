matchit2subclass <- function(treat, X, data, distance, discarded, is.full.mahalanobis,
                             match.matrix=NULL, subclass=6, sub.by="treat",
                             verbose = FALSE){

  if(verbose)
    cat("Subclassifying... \n")  

  # sub.by
  if(!is.vector(sub.by)|length(sub.by)!=1){
    warning(sub.by," is not a valid sub.by option; sub.by=\"treat\" used instead",call.=FALSE); sub.by <- "treat"}
  if(!sub.by%in%c("treat","control","all")){
    warning(sub.by, " is not a valid sub.by option; sub.by=\"treat\" used instead",call.=FALSE); sub.by <- "treat"}

  #subclass
  if(length(subclass)==1){ 
    if(subclass<0 | subclass==1){
      warning(subclass, " is not a valid subclass; subclass=6 used instead",call.=FALSE); subclass <- 6}
  } else {
    if(!is.vector(subclass)){
      warning(subclass, " is not a valid subclass; subclass=6 used instead",call.=FALSE); subclass <- 6}
    if(sum(subclass<=1 & subclass>=0)!=length(subclass)){
      warning("Subclass ", subclass, " is not bounded by 0 and 1; subclass=6 used instead",
              call.=FALSE); subclass <- 6}
  }

  in.sample <- !discarded
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
  p1 <- distance[treat==1][m1]
  p0 <- distance[treat==0][m0]
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
      q <- c(0,quantile(distance,probs=c(subclass)),1)
    else 
      stop("Invalid input for sub.by")
  }
  else {
    if(subclass<=0){stop("Subclass must be a positive vector",call.=FALSE)}
    sprobs <- seq(0,1,length=(round(subclass)+1))
    sprobs <- sprobs[2:(length(sprobs)-1)]
    min.dist <- min(distance,na.rm=TRUE)-0.01 
    max.dist <- max(distance,na.rm=TRUE)+0.01
    if(sub.by=="treat")
      q <- c(min.dist,quantile(p1,probs=sprobs,na.rm=TRUE),
             max.dist)
    else if(sub.by=="control") 
      q <- c(min.dist,quantile(p0,probs=sprobs,na.rm=TRUE),
             max.dist)
    else if(sub.by=="all") 
      q <- c(min.dist,
             quantile(distance,probs=sprobs,na.rm=TRUE),
             max.dist)
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
    psclass <- psclass+i*as.numeric(distance<q2 & distance>=q1)
  }

  ## No subclass for discarded or unmatched units
  psclass[in.sample==0] <- NA
  psclass[!matched] <- NA
  if(verbose){cat("Done\n")}
  res <- list(subclass = psclass, q.cut = q,
              weights = weights.subclass(psclass, treat))

  #warning for discrete data
  unique.classes <- unique(psclass)
  unique.classes <- unique.classes[!is.na (unique.classes)]
  if(length(unique.classes)!=subclass){
    warning("Due to discreteness in data, fewer subclasses generated",call.=F)
  }
 
  class(res) <- c("matchit.subclass", "matchit")
  return(res)
}
