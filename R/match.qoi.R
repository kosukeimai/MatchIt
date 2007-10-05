## Function to calculate summary stats
qoi <- function(xx,tt,ww, t.plot=NULL, c.plot=NULL, sds=NULL,
                standardize = FALSE, std=F){
  weighted.var <- function(x, w) {
    sum(w * (x - weighted.mean(x,w))^2)/(sum(w) - 1)}
  xsum <- matrix(NA,2,7)
  xsum <- as.data.frame(xsum)
  row.names(xsum) <- c("Full","Matched")
  if (standardize)
    names(xsum) <- c("Means Treated","Means Control",
                     "SD Control", "Std. Mean Diff.",
                     "eCDF Med", "eCDF Mean", "eCDF Max")
  else
    names(xsum) <- c("Means Treated","Means Control",
                     "SD Control", "Mean Diff",
                     "eQQ Med", "eQQ Mean", "eQQ Max")
  x1 <- xx[tt==1]
  x0 <- xx[tt==0]
  ww1 <- ww[tt==1]
  ww0 <- ww[tt==0]
  xsum[1,1] <- mean(x1,na.rm=T)
  xsum[1,2] <- mean(x0,na.rm=T)
  xsum[1,3] <- sd(x0,na.rm=T)
  X.t.m <- xx[tt==1][ww1>0]
  X.c.m <- xx[tt==0][ww0>0]
  xsum[2,1] <- weighted.mean(X.t.m, ww1[ww1>0])
  xsum[2,2] <- weighted.mean(X.c.m, ww0[ww0>0])
  xsum[2,3] <- sqrt(weighted.var(X.c.m, ww0[ww0>0]))
  if(!(sum(tt==1)<2|(sum(tt==0)<2))){ 
    xsd <- sd(x1,na.rm=T)
    qqall <- qqsum(x1,x0,standardize=standardize)
    xsum[1,5:7] <- c(qqall$meddiff,qqall$meandiff,qqall$maxdiff)
    if (standardize)
      if (!is.null(sds))
        xsum[1,4] <- (mean(x1,na.rm=T)-mean(x0,na.rm=T))/sds
      else
        xsum[1,4] <- (mean(x1,na.rm=T)-mean(x0,na.rm=T))/xsd
    else
      xsum[1,4] <- mean(x1,na.rm=T)-mean(x0,na.rm=T)
    if(!is.null(t.plot))
      qqmat <- qqsum(xx[t.plot],xx[c.plot],standardize=standardize)
    else
      qqmat <- qqsum(x1[ww1>0],x0[ww0>0],standardize=standardize)
    xsum[2,5:7] <- c(qqmat$meddiff,qqmat$meandiff,qqmat$maxdiff)
    if (standardize)
      if (!is.null(sds))
        xsum[2,4] <- (xsum[2,1]-xsum[2,2])/sds
      else
        xsum[2,4] <- (xsum[2,1]-xsum[2,2])/xsd
    else
      xsum[2,4] <- xsum[2,1]-xsum[2,2]
  } 
  if(!std){
    xsum <- xsum[,c(1:2,4:7)]
  }
  xsum
}

## By subclass
qoi.by.sub <- function(xx,tt,ww,qq,standardize=FALSE){
  qbins <- max(qq,na.rm=TRUE)
  q.table <- matrix(0,6,qbins)
  qn <- matrix(0,3,qbins)
  matched <- ww!=0
  for (i in 1:qbins) {
    qi <- qq[matched]==i & (!is.na(qq[matched]))
    qx <- xx[matched][qi]
    qt <- tt[matched][qi]
    qw <- as.numeric(ww[matched][qi]!=0)
    if(sum(qt==1)<2|(sum(qt==0)<2)){
      if(sum(qt==1)<2)
        warning("Not enough treatment units in subclass ",i,call.=FALSE)
      else if(sum(qt==0)<2)
        warning("Not enough control units in subclass ",i,call.=FALSE)
    }
    qoi.i <- qoi(qx,qt,qw, sds=sd(xx[tt==1],na.rm=T), standardize=standardize)
    q.table[,i] <- as.numeric(qoi.i[1,])
    qn[,i] <- c(sum(qt),sum(qt==0),length(qt))
  }
  q.table <- as.data.frame(q.table)
  qn <- as.data.frame(qn)
  names(q.table) <- names(qn) <- paste("Subclass",1:qbins)
  row.names(q.table) <- names(qoi.i)
  row.names(qn) <- c("Treated","Control","Total")
  list(q.table=q.table,qn=qn)
}
