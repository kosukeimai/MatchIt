matchdef <-  function(formula, in.sample, pscore, nearest = TRUE,
                      replace = FALSE, m.order = 2, ratio = 1,
                      caliper = 0, calclosest = FALSE, mahvars = NULL,
                      exact = FALSE, data = NULL, counter = TRUE, opt
                      = FALSE, ...){ 
  treata <- model.frame(formula,data)[,1,drop=FALSE]
  treat <- as.vector(treata[,1])
  names(treat) <- row.names(treata)
  
  covariates <- model.frame(delete.response(terms(formula)),data)[,,drop=FALSE]
  
  ## Total number of units
  n <- length(treat)
  n1 <- length(treat[treat==1])
  n0 <- length(treat[treat==0])
  
  p1 <- pscore[treat==1]
  p0 <- pscore[treat==0]

  ## Generating separate indices for treated and control units
  if(is.null(names(treat))){names(treat) <- seq(1,n)}
  labels <- names(treat)
  clabels <- labels[treat==0]
  tlabels <- labels[treat==1]
  
  ## Generating match matrix
  match.matrix <- as.data.frame(matrix(0,n1,ratio), row.names = tlabels)

  ## Vectors of whether unit has been matched:
  ## = 0 if not matched (unit # of match if matched)
  ## = -1 if can't be matched (if in.sample=0)
  matchedc <- rep(0,length(p0))
  names(matchedc) <- clabels
  
  ## optimal ratio matching
  if (opt) {
    cat("Optimal Matching Treated: ")
    distance <- matrix(0, ncol=n0, nrow=n1)
    rownames(distance) <- row.names(treata)[treat==1]
    colnames(distance) <- row.names(treata)[treat==0]
    for (i in 1:n1)
      distance[i,] <- abs(p1[i]-p0)
    full <- as.matrix(fullmatch(distance, min.controls = ratio,
                                max.controls = ratio,
                                omit.fraction = (n0-ratio*n1)/n0,...))
    psclass <- full[pmatch(row.names(treata), row.names(full)),]
    psclass <- as.numeric(as.factor(psclass))
    names(psclass) <- row.names(treata)
    for (i in 1:n1) {
      match.matrix[i,] <-
        names(psclass)[match(psclass[tlabels[i]],
                             psclass[-pmatch(tlabels[i],
                                             names(psclass))])]
    }
 
  }
  else if(nearest) {    
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
    sd.cal <- caliper*sqrt(var(pscore[in.sample==1]))
    
    ## Var-covar matrix for Mahalanobis (currently set for full sample)      
    if (!is.null(mahvars)) {
      if(!sum(mahvars%in%names(data))==length(mahvars))
        stop("Mahvars not contained in data",call.=FALSE)
      mahvars <- as.matrix(data[,mahvars])
      row.names(mahvars) <- labels
      Sigma <- var(mahvars)
    }   
    
    ## Now for exact matching
    ## If exact=T, set exact to covariates from distance
    ## Make new variable that is T/F for exact matching, then exact will have the vars if
    ## exactmatch=T
    exactmatch <- exact
    if (!is.logical(exact)) {
      exactmatch <- TRUE} else if (identical(exact,TRUE)) {exact <- covariates}
      
    if (exactmatch==TRUE){
      if(!sum(exact%in%names(data))==length(exact))
        stop("Exact variables not contained in data",call.=FALSE)
      exact <- as.matrix(data[,exact])
      row.names(exact) <- labels
    }   
    ## Looping through nearest neighbour matching for all treatment units
    ## Only do matching for units with in.sample==1 (matched!=-1)
    if(counter==TRUE){
      trseq <- floor(seq(tr/10,tr,tr/10))
      cat("Matching Treated: ")
    }
    
    for(i in 1:tr){
      ## Make new matchedc column to be used for exact matching
      ## Will only be 0 (eligible for matching) if it's an exact match
      if(counter==TRUE) {if(i%in%trseq){cat(10*which(trseq==i),"%...",sep="")}}  # a counter
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
      
      if(m.order==2) {iterp1 <- max(p1[matchedt==0],na.rm=T)}
      if(m.order==3) {iterp1 <- min(p1[matchedt==0],na.rm=T)}
      if(m.order==4) {iterp1 <- sample(p1[matchedt==0][!is.na(p1[matchedt==0])],1)}
      
      ## The treatment unit for this iteration, again resolving ties randomly
      itert <- as.vector(na.omit(tlabels[iterp1==p1 & matchedt==0]))
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
      if (exactmatch==TRUE) {
        for (k in 1:dim(exact)[2]) matchedc2[exact[itert,k]!=exact[clabels,k]] <- -2
      }
      
      ## Need to add a check in case there aren't any eligible matches left...
      if(replace==TRUE & r!=1) {
        if (sum(!clabels%in%match.matrix[itert,(1:r-1)] & matchedc2==0)==0) { 
          deviation <- NULL 
          mindev <- NA
        }
        else
          deviation <- abs(p0[!clabels%in%match.matrix[itert,(1:r-1)] & matchedc2==0]-iterp1)
      }	else { 
        if (sum(matchedc2==0)==0) { 
          deviation <- NULL
          mindev <- NA
        }
        else deviation <- abs(p0[matchedc2==0]-iterp1)
      }
      
      if (caliper!=0 & (!is.null(deviation))) {
        if(replace==TRUE & r!=1)
          pool <- clabels[!clabels%in%match.matrix[itert,(1:r-1)]
                          & matchedc2==0][deviation <= sd.cal]
        else
          pool <- clabels[matchedc2==0][deviation <= sd.cal]
        if(length(pool)==0) { 
          if (calclosest==FALSE) mindev <- NA
          else { 
            if (replace==TRUE & r!= 1){ 
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
      }        else if(!is.null(deviation)) { 
        if (replace==TRUE & r!=1){ 
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
      if (replace==TRUE) matchedc[goodmatch==clabels] <- 0
      
    }
  }
  if(counter==TRUE){cat("Done\n")}
  if(!nearest){match.matrix <- NULL}   else {
    x <- as.matrix(match.matrix)
    x[x==-1] <- NA
    match.matrix <- as.data.frame(x)}
  list(match.matrix = match.matrix, in.sample = in.sample)
}
