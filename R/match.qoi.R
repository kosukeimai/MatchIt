## Function to calculate summary stats
qoi <- function(xx,tt,ww){
  weighted.var <- function(x, w) {
    sum(w * (x - weighted.mean(x,w))^2)/(sum(w) - 1)}
  xsum <- matrix(NA,2,7)
  xsum <- as.data.frame(xsum)
  row.names(xsum) <- c("Full","Matched")
  names(xsum) <- c("Means Treated","Means Control","Pooled SD","QQ Med",
                   "QQ Mean", "QQ Max","Std. Bias")
  x1 <- xx[tt==1]
  x0 <- xx[tt==0]
  ww1 <- ww[tt==1]
  ww0 <- ww[tt==0]
  xsum[1,1] <- mean(x1,na.rm=T)
  xsum[1,2] <- mean(x0,na.rm=T)
  X.t.m <- xx[tt==1][ww1>0]
  X.c.m <- xx[tt==0][ww0>0]
  xsum[2,1] <- weighted.mean(X.t.m, ww1[ww1>0])
  xsum[2,2] <- weighted.mean(X.c.m, ww0[ww0>0])
  if(!(sum(tt==1)<2|(sum(tt==0)<2))){ 
    xsum[1,3] <- sd(xx,na.rm=T)
    qqall <- qqsum(x1,x0)
    xsum[1,4:6] <- c(qqall$meddiff,qqall$meandiff,qqall$maxdiff)
    xsum[1,7] <- (mean(x1,na.rm=T)-mean(x0,na.rm=T))/sd(x1,na.rm=T)
    xsum[2,3] <- sqrt(weighted.var(xx[ww>0],ww[ww>0]))
    qqmat <- qqsum(x1[ww1>0],x0[ww0>0])
    xsum[2,4:6] <- c(qqmat$meddiff,qqmat$meandiff,qqmat$maxdiff)
    xsum[2,7] <- (xsum[2,1]-xsum[2,2])/sqrt(weighted.var(X.t.m,ww1[ww1>0]))
  } 
  xsum
}
