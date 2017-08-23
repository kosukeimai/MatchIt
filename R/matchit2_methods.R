# matchit2cem - matchit wrapper for cem matching algorithm
#
# 06/10/2008 - m.blackwell
#
# this function takes inputs from matchit() and returns the
# strata for each observation in the subclass entry and the
# weight for each observation in the weight entry. No match
# matrix is returned since matches are not unique within
# strata. 
#
#----------------------------------------------------------
# MATCHIT method= cem
#----------------------------------------------------------
matchit2cem <- function(treat, X, data, distance, discarded, is.full.mahalanobis,
                            ratio = 1, verbose = FALSE, k2k.method=NULL, ...) {

  
  if (!requireNamespace("cem", quietly = TRUE)) 
      stop("cem package is required.  Please install it.")

  if (verbose)
    cat("Coarsened exact matching...\n")
  
  n <- length(treat)
  
  # cem takes the data all together and wants the treatment specified
  # with the column name of the data frame. Here we massage the matchit
  # inputs to this format. Note that X has its proper columnames, but
  # treat does not have the original column name. 
  cem.data <- as.data.frame(cbind(treat,X))
  
  mat <-
    cem::cem(treatment="treat",data=cem.data,verbose=as.integer(verbose)+1,
        method=k2k.method,...)

  # here we create a column vector where the matched entry get its stratum
  # and the unmatched entry gets an NA.
  strat <- rep(NA,n)
  names(strat) <- names(treat)
  strat[mat$matched] <- mat$strata[mat$matched]

  # here we just add the names onto the wieght from the cem output
  wh <- mat$w
  names(wh) <- names(treat)

  # weighting functions in matchit error-out on these conditions,
  # so we should too.
 
  if (sum(wh)==0) 
    stop("No units were matched")
  else if (sum(wh[treat==1])==0)
    stop("No treated units were matched")
  else if (sum(wh[treat==0])==0)
    stop("No control units were matched")
  
  res <- list(subclass = strat, weights = mat$w)
  class(res) <- "matchit"
  return(res)
}

#----------------------------------------------------------
# MATCHIT method= exact
#----------------------------------------------------------
matchit2exact <- function(treat, X, data, distance, discarded, is.full.mahalanobis, verbose=FALSE, ...){
  
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

#----------------------------------------------------------
# MATCHIT method= full
#----------------------------------------------------------
matchit2full <- function(treat, X, data, distance, discarded, is.full.mahalanobis,
                         verbose=FALSE, ...) { 
  
  if(verbose)
    cat("Full matching... \n")
  
  ## full matching for undiscarded units
  ttt <- treat[!discarded]
  ddd <- distance[!discarded]
  n0 <- length(ttt[ttt==0])
  n1 <- length(ttt[ttt==1])
  
  if (is.matrix(distance)){
    if (ncol(distance) != length(treat) | nrow(distance) != length(treat))
      error("dimension of distance matrix is incorrect")
    d <- distance
    d <- d[treat == 1 & !discarded, treat == 0 & !discarded]
  } else {
    d1 <- ddd[ttt==1]
    d0 <- ddd[ttt==0]
    d <- matrix(0, ncol=n0, nrow=n1)
    rownames(d) <- names(ttt[ttt==1])
    colnames(d) <- names(ttt[ttt==0])
    for (i in 1:n1) 
        d[i,] <- abs(d1[i]-d0)
  }     
  
  full <- optmatch::fullmatch(d, ...)
  psclass <- full[pmatch(names(ttt), names(full))]
  psclass <- as.numeric(as.factor(psclass))
  names(psclass) <- names(ttt)
  
  ## add psclass = NA for discarded units
  if (any(discarded)) {
    tmp <- rep(NA, sum(discarded))
    names(tmp) <- names(treat[discarded])
    psclass <- c(psclass, tmp)[names(treat)]
  }
  
  ## calculate weights and return the results
  res <- list(subclass = psclass, weights = weights.subclass(psclass, treat))
  class(res) <- c("matchit.full", "matchit")
  return(res)
}

#----------------------------------------------------------
# MATCHIT method= genetic
#----------------------------------------------------------
matchit2genetic <- function(treat, X, data, distance, discarded, is.full.mahalanobis,
                            ratio = 1, verbose = FALSE, ...) {
  
  if (verbose)
    cat("Genetic matching... \n")
  
  tt <- treat[!discarded]
  n <- length(tt)
  n1 <- length(tt[tt==1])
  xx <- X[!discarded,]
  dd <- distance[!discarded]
  tind <- (1:n)[tt==1]
  cind <- (1:n)[tt==0]
  labels <- names(tt)
  tlabels <- names(tt[tt==1])
  clabels <- names(tt[tt==0])
  out <- Matching::GenMatch(tt, cbind(dd, xx), M = ratio, ...)$matches
  ## ratio matching does not seem to work with GenMatch
  mm <- matrix(0, nrow = n1, ncol = max(table(out[,1])), dimnames =
                 list(tlabels, 1:max(table(out[,1]))))
  for (i in 1:n1) {
    tmp <- labels[c(out[out[,1]==tind[i],2:(ratio+1)])]
    if (length(tmp) < ncol(mm))
      tmp <- c(tmp, rep(NA, ncol(mm)-length(tmp)))
    mm[i,] <- tmp
  }
  
  if (any(discarded)) {
    tdisc <- discarded[treat==1]
    tmp <- matrix(NA, nrow = sum(tdisc), ncol = ncol(mm), dimnames =
                    list(names(treat[treat == 1 & discarded]),
                         1:ncol(mm)))
    mm <- as.matrix(rbind(mm, tmp)[names(treat[treat==1]),])
  }
  
  res <- list(match.matrix = mm, weights = weights.matrix(mm, treat,
                                                          discarded))
  class(res) <- "matchit"
  return(res)
}

#----------------------------------------------------------
# MATCHIT method= nearest
#----------------------------------------------------------
matchit2nearest <-  function(treat, X, data, distance, discarded,
                             ratio=1, replace = FALSE, m.order = "largest",  
                             caliper = 0, calclosest = FALSE,
                             mahvars = NULL, exact = NULL,
                             subclass=NULL, verbose=FALSE, sub.by=NULL,
                             is.full.mahalanobis,...){  
  
  if(verbose)
    cat("Nearest neighbor matching... \n")
  
  #replace
  if(!(identical(replace,TRUE) | identical(replace,FALSE))){
    warning("replace=",replace," is invalid; used replace=FALSE instead",call.=FALSE);replace=FALSE}
  #m.order
  if(!(identical(m.order,"largest") | identical(m.order,"smallest") |
       identical(m.order,"random"))){
    warning("m.order=",m.order," is invalid; used m.order='largest' instead",call.=FALSE);m.order="largest"}
  #ratio
  ratio <- round(ratio)
  if(!is.numeric(ratio) | ratio[1]<1 | !identical(round(length(ratio)),1)){
    warning("ratio=",ratio," is invalid; used ratio=1 instead",call.=FALSE);ratio=1}
  #caliper
  if(!is.vector(caliper) | !identical(round(length(caliper)),1)){
    warning("caliper=",caliper," is invalid; Caliper matching not done",call.=FALSE);caliper=0}
  if(caliper<0){
    warning("caliper=",caliper," is less than 0; Caliper matching not done",call.=FALSE);caliper=0}
  #calclosest
  if(!(identical(calclosest,TRUE)| identical(calclosest,FALSE))){
    warning("calclosest=",calclosest," is invalid; used calclosest=FALSE instead",call.=FALSE)
    calclosest=FALSE}
  #mahvars & caliper
  if (!is.null(mahvars) & caliper[1]==0){
    warning("No caliper size specified for Mahalanobis matching.  Caliper=.25 used.",call. = FALSE);caliper=.25}
  #when mahalanobis distance is used for all covars
  if(is.full.mahalanobis){
    mahvars <- X
    Sigma <- var(X)
    ## Note: caliper irrelevant, but triggers mahalanobis matching
    caliper <- .25
    ## no subclass with full mahalanobis
    if(!is.null(subclass)){
      warning("No subclassification with pure Mahalanobis distance.",call. = FALSE)
      subclass <- NULL
    }
  }
  
  # Sample sizes, labels
  n <- length(treat)
  n0 <- length(treat[treat==0])
  n1 <- length(treat[treat==1])
  d1 <- distance[treat==1]
  d0 <- distance[treat==0]
  
  if(is.null(names(treat)))
    names(treat) <- 1:n
  labels <- names(treat)
  tlabels <- names(treat[treat==1])
  clabels <- names(treat[treat==0])
  in.sample <- !discarded
  names(in.sample) <- labels
  
  ## 10/1/07: Warning for if fewer control than ratio*treated and matching without replacement
  if (n0 < ratio*n1 & replace==FALSE) {
    if (ratio > 1)  warning(paste("Not enough control units for ", ratio, " matches for each treated unit when matching without replacement.  Not all treated units will receive", ratio, "matches"))
    else warning(paste("Fewer control than treated units and matching without replacement.  Not all treated units will receive a match.  Treated units will be matched in the order specified by m.order:", m.order))
  }
  
  ## Generating match matrix
  match.matrix <- matrix(0, nrow=n1, ncol=ratio, dimnames=list(tlabels, 1:ratio))
  
  ## Vectors of whether unit has been matched:
  ## = 0 if not matched (unit # of match if matched)
  ## = -1 if can't be matched (if in.sample=0)
  matchedc <- rep(0,length(d0))
  names(matchedc) <- clabels
  
  ## These are the units that are ineligible because of discard
  ## (in.sample==0) 
  matchedc[in.sample[clabels]==0] <- -1
  match.matrix[in.sample[tlabels]==0,] <- -1
  matchedt <- match.matrix[,1] 
  names(matchedt) <- tlabels
  
  ## total number of matches (including ratios) = ratio * n1
  tr <- length(match.matrix[match.matrix!=-1])
  r <- 1
  
  ## Caliper for matching (=0 if caliper matching not done)
  sd.cal <- caliper*sqrt(var(distance[in.sample==1]))
  
  ## Var-covar matrix for Mahalanobis (currently set for full sample)
  if (!is.null(mahvars) & !is.full.mahalanobis) {
    if(!sum(mahvars%in%names(data))==length(mahvars)) {
      warning("Mahvars not contained in data.  Mahalanobis matching not done.",call.=FALSE)
      mahvars=NULL
    }
    else {  ww <- mahvars%in%dimnames(X)[[2]]
    nw <- length(mahvars)
    mahvars <- data[,mahvars,drop=F]
    Sigma <- var(mahvars)
    if(sum(ww)!=nw){
      X <- cbind(X,mahvars[!ww])
    }
    mahvars <- as.matrix(mahvars)
    }
  }
  
  ## Now for exact matching within nearest neighbor
  ## exact should not equal T for this type of matching--that would get sent to matchit2exact
  if (!is.null(exact)){
    if(!sum(exact%in%names(data))==length(exact)) {
      warning("Exact variables not contained in data. Exact matching not done.",call.=FALSE)
      exact=NULL
    }
    else {
      ww <- exact%in%dimnames(X)[[2]]
      nw <- length(exact)
      exact <- data[,exact,drop=F]
      if(sum(ww)!=nw){
        X <- cbind(X,exact[!ww])
      }
    }
  }
  
  ## Looping through nearest neighbour matching for all treatment units
  ## Only do matching for units with in.sample==1 (matched!=-1)
  if(verbose){
    trseq <- floor(seq(tr/10,tr,tr/10))
    cat("Matching Treated: ")
  }
  
  for(i in 1:tr){
    ## Make new matchedc column to be used for exact matching
    ## Will only be 0 (eligible for matching) if it's an exact match
    if(verbose) {if(i%in%trseq){cat(10*which(trseq==i),"%...",sep="")}}  # a counter
    matchedc2 <- matchedc
    ##in cases there's no replacement and all controls have been used up
    if(!0%in%matchedc2){  
      match.matrix[match.matrix[,r]==0 & !is.na(match.matrix[,r]),r] <- NA
      if(r<ratio){match.matrix[,(r+1):ratio] <- NA}
      break
    }
    
    ##in case there's replacement, but all units have been used in
    ##previous ratios
    if(sum(!is.na(match.matrix[,r]))==0){
      if(r<ratio){match.matrix[,(r+1):ratio] <- NA}
      break
    }
    
    ## Which ratio we're on
    if(r!=ceiling(i/(tr/ratio))) {r <- r+1; matchedt <- match.matrix[,r]}
    
    if(m.order=="largest") {iterd1 <- max(d1[matchedt==0],na.rm=T)}
    if(m.order=="smallest") {iterd1 <- min(d1[matchedt==0],na.rm=T)}
    if(m.order=="random") {iterd1 <- sample(d1[matchedt==0][!is.na(d1[matchedt==0])],1)}
    
    ## The treatment unit for this iteration, again resolving ties randomly
    itert <- as.vector(na.omit(tlabels[iterd1==d1 & matchedt==0]))
    if(length(itert)>1){itert <- sample(itert,1)}
    
    ## Calculating all the absolute deviations in propensity scores
    ## Calculate only for those eligible to be matched (matchedc==0)
    ## this first if statement only applies to replacement ratio
    ## matching, so that each treatment unit is matched to a different
    ## control unit than from the previous round
    
    ## match number = NA if no units within caliper
    
    ## Set things up for exact matching
    ## Make matchedc2==-2 if it isn't an exact match
    ## There might be a more efficient way to do this, but I couldn't figure
    ## out another way to compare a vector with the matrix
    if (!is.null(exact)) {
      for (k in 1:dim(exact)[2]) matchedc2[exact[[k]][itert]!=exact[[k]][clabels]] <- -2
    }
    
    ## Need to add a check in case there aren't any eligible matches left...
    if(replace & r!=1) {
      if (sum(!clabels%in%match.matrix[itert,(1:r-1)] & matchedc2==0)==0) { 
        deviation <- NULL 
        mindev <- NA
      }
      else
        deviation <- abs(d0[!clabels%in%match.matrix[itert,(1:r-1)] & matchedc2==0]-iterd1)
    }
    else { 
      if (sum(matchedc2==0)==0) { 
        deviation <- NULL
        mindev <- NA
      }
      else deviation <- abs(d0[matchedc2==0]-iterd1)
    }
    
    if (caliper!=0 & (!is.null(deviation))) {
      if(replace & r!=1)
        pool <- clabels[!clabels%in%match.matrix[itert,(1:r-1)]
                        & matchedc2==0][deviation <= sd.cal]
      else
        pool <- clabels[matchedc2==0][deviation <= sd.cal]
      if(length(pool)==0) { 
        if (calclosest==FALSE) mindev <- NA
        else { 
          if (replace & r!= 1){ 
            mindev <- clabels[!clabels%in%match.matrix[itert,(1:r-1)]][min(deviation)==deviation]
          } else{mindev <- clabels[matchedc2==0][min(deviation)==deviation]}
        }
      }
      else if (length(pool)==1) mindev <- pool[1]
      else if (is.null(mahvars)) mindev <- sample(pool, 1)
      else {
        ## This has the important vars for the C's within the caliper
        poolvarsC <- mahvars[pool,,drop=F]
        ## Sigma is the full group var/covar matrix of Mahalvars
        mahal <- mahalanobis(poolvarsC, mahvars[itert,],Sigma)
        mindev <- pool[mahal==min(mahal)]
      }
    }
    else if(!is.null(deviation)) { 
      if (replace & r!=1){ 
        mindev <- clabels[!clabels%in%match.matrix[itert,(1:r-1)] & matchedc2==0][min(deviation)==deviation]
      } else {mindev <- clabels[matchedc2==0][min(deviation)==deviation]}
    }
    
    ## Resolving ties in minimum deviation by random draw
    if(length(mindev)>1){goodmatch <- sample(mindev,1)} else goodmatch <- mindev
    
    ## Storing which treatment unit has been matched to control, and
    ## vice versa
    matchedt[itert==tlabels] <- goodmatch
    matchedc[goodmatch==clabels] <- itert
    
    ## instead of the in.sample, we now have an index with dimensions n1 by # of
    ## matches (ratio)
    match.matrix[which(itert==tlabels),r] <- goodmatch
    
    ## If matching with replacement, set matchedc back to 0 so it can be reused
    if (replace) matchedc[goodmatch==clabels] <- 0
    
  }
  if(verbose){cat("Done\n")}
  
  x <- as.matrix(match.matrix)
  x[x==-1] <- NA
  
  ## Calculate weights and return the results
  res <- list(match.matrix = match.matrix, weights =
                weights.matrix(match.matrix, treat, discarded), X=X)
  
  ## Subclassifying
  if(!is.null(subclass)){
    if(is.null(sub.by)) sub.by="treat"
    psres <- matchit2subclass(treat,X,data,distance,discarded,
                              match.matrix=match.matrix,
                              subclass=subclass,
                              verbose=verbose, sub.by=sub.by, ...)
    res$subclass <- psres$subclass
    res$q.cut <- psres$q.cut
    class(res) <- c("matchit.subclass", "matchit")
  } else{
    class(res) <- "matchit"
  }
  return(res)
}

#----------------------------------------------------------
# MATCHIT method= optimal
#----------------------------------------------------------
matchit2optimal <- function(treat, X, data, distance, discarded, is.full.mahalanobis, 
                            ratio = 1, verbose=FALSE, ...) {
  
  if (!requireNamespace("optmatch", quietly = TRUE)) 
    stop("optmatch package is required.  Please install it.")
  
  if(verbose)
    cat("Optimal matching... \n")
  
  ## optimal matching for undiscarded units
  ttt <- treat[!discarded]
  n0 <- length(ttt[ttt==0])
  n1 <- length(ttt[ttt==1])
  d1 <- distance[ttt==1]
  d0 <- distance[ttt==0]
  d <- matrix(0, ncol=n0, nrow=n1)
  tlabels <- rownames(d) <- names(ttt[ttt==1])
  clabels <- colnames(d) <- names(ttt[ttt==0])
  for (i in 1:n1) 
    d[i,] <- abs(d1[i]-d0)
  full <- optmatch::fullmatch(d, min.controls = ratio,
                              max.controls = ratio,
                              omit.fraction = (n0-ratio*n1)/n0, ...)
  psclass <- full[pmatch(names(ttt), names(full))]
  psclass <- as.numeric(as.factor(psclass))
  names(psclass) <- names(ttt)
  
  mm <- matrix(0, nrow = n1, ncol = ratio, dimnames = list(tlabels, 1:ratio))
  for (i in 1:n1)
    mm[i,] <- names(which(psclass[tlabels[i]] == psclass[-pmatch(tlabels[i],
                                                                 names(psclass))]))
  
  if (any(discarded)) {
    ## add psclass = NA for discarded units
    tmp <- rep(NA, sum(discarded))
    names(tmp) <- names(treat[discarded])
    psclass <- c(psclass, tmp)[names(treat)]
    
    ## add match.matrix = NA for discarded units
    tdisc <- discarded[treat==1]
    if (any(tdisc)) {
      tmp <- matrix(NA, nrow = sum(tdisc), ncol= ratio,
                    dimnames = list(names(treat[treat==1 & discarded]),
                                    1:ratio))
      mm <- as.matrix(rbind(mm, tmp)[names(treat[treat==1]),])
    }
  }
  
  ## calculate weights and return the results
  res <- list(match.matrix = mm, subclass = psclass,
              weights = weights.matrix(mm, treat, discarded))
  
  class(res) <- "matchit"
  return(res)
}

#----------------------------------------------------------
# MATCHIT method= subclass
#----------------------------------------------------------
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
