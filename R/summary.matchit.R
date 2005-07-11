summary.matchit <- function(obj, verbose=F, ...) {
  weighted.var <- function(x, w) {
    sum(w * (x - weighted.mean(x,w))^2)/(sum(w) - 1)}
  ## Function to calculate summary stats
  qoi <- function(xx,tt,ww){
    xsum <- matrix(0,2,7)
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
    if(sum(tt==1)<2|(sum(tt==0)<2)){
      xsum[c(1,2),c(3,4,5)] <- NA
    } else {
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
  XX <- foo$X
  treat <- foo$treat
  weights <- foo$weights
  aa <- apply(XX,2,qoi,tt=treat,ww=weights)
  nam <- dimnames(obj$X)[[2]]
  kk <- ncol(obj$X)
  sum.all <- as.data.frame(matrix(0,kk,7))
  sum.matched <- as.data.frame(matrix(0,kk,7))
  row.names(sum.all) <- row.names(sum.matched) <- nam
  names(sum.all) <- names(sum.matched) <- names(aa[[1]])
  sum.all.int <- sum.matched.int <- NULL
  for(i in 1:kk){
    sum.all[i,] <- aa[[i]][1,]
    sum.matched[i,] <- aa[[i]][2,]
    if(verbose){
      for(j in i:kk){
        x2 <- XX[,i]*as.matrix(XX[,j])
        jqoi <- qoi(x2,tt=treat,ww=weights)
        sum.all.int <- rbind(sum.all.int,jqoi[1,])
        sum.matched.int <- rbind(sum.matched.int,jqoi[2,])
        row.names(sum.all.int)[nrow(sum.all.int)] <-
          row.names(sum.matched.int)[nrow(sum.matched.int)] <-
            paste(nam[i],nam[j],sep="x")
      }
    }
  }
  xn <- aa[[1]]$xn
  sum.all <- rbind(sum.all,sum.all.int)
  sum.matched <- rbind(sum.matched,sum.matched.int)
  nn <- rbind(table(obj$treat),
              table(obj$weights!=0,obj$treat)[2:1,])
  dimnames(nn) <- list(c("Full","Matched","Discarded"),
                       c("Control","Treated"))
  obj$nn <- nn
  obj$sum.all <- sum.all
  obj$sum.matched <- sum.matched
  obj$verbose <- verbose
  class(obj) <- "summary.matchit"
  return(obj)
}  
