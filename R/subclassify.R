subclassify <- function(formula,data,in.sample,pscore,nearest=TRUE,
                        match.matrix,subclass=0,sub.by="treat",
                        counter=TRUE, full = FALSE, full.options=list()){

  data <- eval(data,parent.frame())
  treata <- model.frame(formula,data)[,1,drop=FALSE]
  treat <- as.vector(treata[,1])
  names(treat) <- row.names(treata)

  if(full) { # full matching with propensity score
    if(counter){cat("Full Matching...")}
    full <- full.options
    if(!is.list(full.options)){
      warning("full.options must be a list; assuming defaults for full matching",call.=FALSE)
    } 
    if(is.null(full$min.controls)){
      full$min.controls <- 0
    }
    if(is.null(full$max.controls)){
      full$max.controls <- Inf
    }
    if(is.null(full$omit.fraction)){
      full$omit.fraction <- NULL
    }
    if(is.null(full$tol)){
      full$tol <- 0.01
    }
    if(is.null(full$subclass.indices)){
      full$subclass.indices <- NULL
    }
    notin <- names(full.options)[which(!names(full.options)%in%c("min.controls","max.controls",
                                "omit.fraction","omit.fraction",
                                "tol", "subclass.indices"))]
    if(!is.null(notin)){
      warning(paste(notin,collapse=" "), " in full.options invalid and ignored for full matching",call.=FALSE)
    }
    n1 <- length(treat[treat==1])
    n0 <- length(treat[treat==0])
    p1 <- pscore[treat==1]
    p0 <- pscore[treat==0]
    distance <- matrix(0, ncol=n0, nrow=n1)
    rownames(distance) <- row.names(treata)[treat==1]
    colnames(distance) <- row.names(treata)[treat==0]
    for (i in 1:n1)
      distance[i,] <- abs(p1[i]-p0)
    full <- as.matrix(fullmatch(distance,subclass.indices = full$subclass.indices,
                                min.controls = full$min.controls,
                                max.controls = full$max.controls,
                                omit.fraction = full$omit.fraction,
                                tol = full$tol))
    psclass <- full[pmatch(row.names(treata), row.names(full)),]
    psclass <- as.numeric(as.factor(psclass))
    names(psclass) <- row.names(treata)
    q <- NULL
    if(counter){cat("Done\n")}
  }
  else if(subclass) {
    if(counter){cat("Subclassifying...")}  
    n <- length(treat)
    if(nearest){
      match.matrix <- match.matrix[in.sample[treat==1],,drop=F]
      t.units <- row.names(match.matrix)[!is.na(match.matrix)]
      
      c.units <- na.omit(as.vector(as.matrix(match.matrix)))
      matched <-c(t.units,c.units)
      matched <- names(treat)%in%matched
    } else{
      matched <- rep(TRUE,n)}
    names(matched) <- names(treat)
    m1 <- matched[treat==1]
    m0 <- matched[treat==0]
    p1 <- pscore[treat==1][m1]
    p0 <- pscore[treat==0][m0]
    
    if(length(subclass)!=1 | (length(subclass)==1 &
               identical(subclass<1,TRUE))) {
      subclass <- sort(subclass)
      if (subclass[1]==0) subclass <- subclass[-1]
      if (subclass[length(subclass)]==1) subclass <- subclass[-length(subclass)]
      if(sub.by=="treat") {
        q <- c(0,quantile(p1,probs=c(subclass)),1)
      }
      else if(sub.by=="control") {
        q <- c(0,quantile(p0,probs=c(subclass)),1)
      }
      else if(sub.by=="all") {
        q <- c(0,quantile(pscore,probs=c(subclass)),1)
      }
      ##          else {stop("Must specify a valid sub.by",call.=FALSE)
      ##             }
    }
    else {
      ##        if(subclass<=0){stop("Subclass must be a positive vector",call.=FALSE)}
      sprobs <- seq(0,1,length=(round(subclass)+1))
      sprobs <- sprobs[2:(length(sprobs)-1)]
      if(sub.by=="treat")            {
        q <- c(0,quantile(p1,probs=sprobs,na.rm=TRUE),1)
      }
      else if(sub.by=="control") {
        q <- c(0,quantile(p0,probs=sprobs,na.rm=TRUE),1)
      }
      else if(sub.by=="all") {
        q <- c(0,quantile(pscore,probs=sprobs,na.rm=TRUE),1)
      }
      ##      else {stop("Must specify a valid sub.by",call.=FALSE)}
    }
    
    qbins <- length(q)-1
    psclass <- rep(0,n)
    names(psclass) <- names(treat)
    for (i in 1:qbins){
      q1 <- q[i]
      q2 <- q[i+1]
      psclass <- psclass+i*as.numeric(pscore<q2 & pscore>=q1)}
    
    ## Make sure not to assign subclass to discarded units
    psclass[in.sample==0] <- 0
    psclass[!matched] <- 0
    if(counter){cat("Done\n")}
  }  
  else {psclass <- NULL; q=NULL}

  return(list(psclass = psclass, q.cut = q))
}
