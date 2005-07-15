neyman.matchit<-function(Y, object, bootstrap=NULL, counter=TRUE){

  if(is.null(bootstrap)) {
    res<-neyman.est(Y, object)
  } else {
    ate.est<-Ntrt<-Ncont<-NULL
    cl<-object$call
    cl$counter<-FALSE
    cl$data <- as.name("dta")
    if(counter) {
      bseq<-floor(seq(bootstrap/10, bootstrap, bootstrap/10))
    }
    for(i in 1:bootstrap){
      sidx <- sample(c(1:nrow(data)), size=nrow(data),
                     replace=T)
      dta <- eval(object$call$data, sys.parent())[sidx,]
      bobject<-eval(cl)
      YY <- Y[sidx]
      tmp<-neyman.est(YY,bobject)
      ate.est<-rbind(ate.est, tmp$ate.est)
      Ntrt<-rbind(Ntrt, tmp$Ntrt)
      Ncont<-rbind(Ncont, tmp$Ncont)
      if(counter & i%in%bseq)
        cat(10*which(bseq==i),"%...", sep="")
    }
    cat("Done\n")
    row.names(ate.est)<-c(1:bootstrap)
    res<-list(ate.est=ate.est, ate.var=NULL, Ntrt=Ntrt, Ncont=Ncont, bootstrap=bootstrap)
  }
    
  # Subclass estimates calculated above...now aggregate
  if(is.null(object$call$sub.by)){
    sub.by <- "treat"
  } else {
    sub.by <- object$call$sub.by
  }
  if(!is.null(bootstrap)){
    if(is.null(object$call$subclass) | identical(object$call$subclass,0)){
      ate<-apply(res$ate.est, 2, mean)
      numtrt <- apply(res$Ntrt, 2, mean)
      numcont <- apply(res$Ncont, 2, mean)
    } else{
      if(sub.by=="treat"){
        w<-res$Ntrt/apply(res$Ntrt,1,sum)
      } else if(sub.by=="control"){
        w<-res$Ncont/apply(res$Ncont,1,sum)
      } else if(sub.by=="all"){
        w <- (res$Ncont+res$Ntrt)/apply((res$Ncont+res$Ntrt),1,sum)
      }
                                        # Set w to 0 if block has fewer then 2 t or c units
      w[res$Ntrt<=1 | res$Ncont<=1] <- rep(0, sum(res$Ntrt<=1 | res$Ncont<=1)) 
      ate.boot <- NULL
      for(i in 1:bootstrap){
        ate.boot<-c(ate.boot,weighted.mean(res$ate.est[i,],w[i,]))
      }
      ate<-c(mean(ate.boot),apply(res$ate.est, 2, mean))
      numtrt <- c(mean(apply(res$Ntrt, 1, sum)), apply(res$Ntrt, 2, mean))
      numcont <- c(mean(apply(res$Ncont, 1, sum)), apply(res$Ncont, 2, mean))
    }
  } else {
    if(sub.by=="treat"){
      w<-res$Ntrt/sum(res$Ntrt)
    } else if(sub.by=="control"){
      w<-res$Ncont/sum(res$Ncont)
    } else if(sub.by=="all"){
      w <- (res$Ncont+res$Ntrt)/sum(res$Ncont+res$Ntrt)
    } 
                                        # Set w to 0 if block has fewer then 2 t or c units
    w[res$Ntrt<=1 | res$Ncont<=1] <- rep(0, sum(res$Ntrt<=1 | res$Ncont<=1)) 
    ate<-apply(res$ate.est, 2, weighted.mean, w)
    if (length(res$Ntrt)> 1) {
      numtrt <- c(sum(res$Ntrt), res$Ntrt)
      numcont <- c(sum(res$Ncont), res$Ncont)
    } else {
      numtrt <- res$Ntrt
      numcont <- res$Ncont
    }
  }
  if(is.null(res$ate.var)){
    if(is.null(object$call$subclass) | identical(object$call$subclass,0)){
      se<-sqrt(apply(res$ate.est, 2, var))
    } else {
      se <- c(sd(ate.boot),sqrt(apply(res$ate.est, 2, var)))
    }
  } else{
    se<-sqrt(res$ate.var)
  }
  if(is.null(bootstrap)){
    if(length(ate)>1){
      if(is.null(res$ate.var)){
        tmp<-apply(res$ate.est, 1, weighted.mean, w)
        ate<-c(mean(tmp), ate)
        se<-c(sqrt(var(tmp)), se)
      }    
      else{
	  # weight=0 if block has < 2 units in either group
        ate<-c(weighted.mean(ate, w, na.rm=T), ate)
	  # Only calculate se for blocks with > 1 unit (se!=nan)	
        se<-c(sqrt(sum((w[!is.nan(se)]*se[!is.nan(se)])^2)), se)
      }
    }
  }

  ressum<-list(ate=ate, se=se, Ntrt=numtrt, Ncont=numcont)
  class(ressum) <- "neyman"
  ressum
}
