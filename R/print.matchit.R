print.matchit <- function(x,  digits = max(3, getOption("digits") -
  3), ...)
  {
    t.test.wtd <- function(x, treat, weights) {
      trt <- cov.wt(as.matrix(x[treat==1]), wt=weights[treat==1])
      mean.t <- trt$center
      cov.t <- trt$cov    
      cont <- cov.wt(as.matrix(x[treat==0]), wt=weights[treat==0])
      mean.c <- cont$center
      cov.c <- cont$cov
      n.t <- sum(weights[treat==1])
      n.c <- sum(weights[treat==0])
      ratio.t <- cov.t/n.t
      ratio.c <- cov.c/n.c
      tt  <- (mean.t-mean.c)/sqrt(ratio.t+ratio.c)
    }
    weighted.var <- function(x, w) {
      sum(w * (x - weighted.mean(x,w))^2)/(sum(w) - 1)}
    qoi <- function(xx,tt,ww){
      xsum <- matrix(0,2,5)
      xsum <- as.data.frame(xsum)
      row.names(xsum) <- c("Full","Matched")
      names(xsum) <- c("Means Treated","Means Control","SD","T-stat","Bias")
      xn <- matrix(0,2,3)
      xn <- as.data.frame(xn)
      row.names(xn) <- c("Full","Matched")
      names(xn) <- c("Treated","Control","Total")
      xn[1,1] <- sum(tt==1)
      xn[1,2] <- sum(tt==0)
      xn[1,3] <- length(tt)
      x1 <- xx[tt==1]
      x0 <- xx[tt==0]
      xsum[1,1] <- mean(x1,na.rm=T)
      xsum[1,2] <- mean(x0,na.rm=T)
      X.t.m <- xx[tt==1][ww[tt==1]>0]
      X.c.m <- xx[tt==0][ww[tt==0]>0]
      xsum[2,1] <- weighted.mean(X.t.m, ww[tt==1][ww[tt==1]>0])
      xsum[2,2] <- weighted.mean(X.c.m, ww[tt==0][ww[tt==0]>0])
      if(sum(tt==1)<2|(sum(tt==0)<2)){
        xsum[c(1,2),c(3,4,5)] <- NA
      } else {
        xsum[1,3] <- sd(xx,na.rm=T) 
        xsum[1,4] <- -1*t.test(xx~tt)$sta
        xsum[1,5] <- (mean(x1,na.rm=T)-mean(x0,na.rm=T))/sd(x1,na.rm=T)
        xsum[2,3] <- sqrt(weighted.var(xx[ww>0],ww[ww>0]))
        xsum[2,4] <- t.test.wtd(xx[ww>0],tt[ww>0],ww[ww>0])
        xsum[2,5] <- (xsum[2,1]-xsum[2,2])/sqrt(weighted.var(X.t.m,ww[tt==1][ww[tt==1]>0]))
      } 
      xn[2,1] <- sum(tt[ww>0]==1)
      xn[2,2] <- sum(tt[ww>0]==0)
      xn[2,3] <- length(tt[ww>0])
      list(xsum=xsum,xn=xn)
    }
    qoi.by.sub <- function(xx,tt,ww,qq){
      qbins <- max(qq,na.rm=TRUE)
      q.table <- matrix(0,5,qbins)
      qn <- matrix(0,3,qbins)
      matched <- ww!=0
      for (i in 1:qbins)
        {
          qi <- qq[matched]==i
          qx <- xx[matched][qi]
          qt <- tt[matched][qi]
          qw <- ww[matched][qi]
          if(sum(qt==1)<2|(sum(qt==0)<2)){
            if(sum(qt==1)<2){warning("Not enough treatment units in subclass ",i,call.=FALSE)}
            else if(sum(qt==0)<2){warning("Not enough control units in subclass ",i,call.=FALSE)}
          }
          qoi.i <- qoi(qx,qt,qw)
          q.table[,i] <- as.numeric(qoi.i$xsum[1,])
          qn[,i] <- as.numeric(qoi.i$xn[1,])
        }
      q.table <- as.data.frame(q.table)
      qn <- as.data.frame(qn)
      names(q.table) <- names(qn) <- paste("Subclass",1:qbins)
      row.names(q.table) <- names(qoi.i$xsum)
      row.names(qn) <- c("Treated","Control","Total")
      list(q.table=q.table,qn=qn)
    }
    
    pscore <- x$data[,"pscore"]
    treat <- eval(x$treat,x$data)
    names(treat) <- names(pscore)
    cat("\nAssignment model specification:\n", deparse(x$call),"\n\n",sep = "")
    if(identical(eval(x$call$exact),TRUE)) #for exact matching
      {
        cat("\nSample sizes for full and exactly matched data:\n\n")
        xqoi <- qoi(pscore, treat, x$psweights)
        print.data.frame(xqoi$xn, digits = digits)
        cat("\n\n")
      } else {        
        if(is.null(x$match.matrix))
          {
            cat("Summary of propensity score for full sample:\n\n")
            xqoi <- qoi(pscore, treat, x$psweights)
            print.data.frame(xqoi$xsum[1,], digits = digits)
            cat("\nSample sizes:\n\n")
            print.data.frame(xqoi$xn[1,], digits = digits)
          } else {
            cat("\nSummary of propensity score for full and matched samples:\n\n")
            xqoi <- qoi(pscore, treat, x$psweights)
            print.data.frame(xqoi$xsum, digits = digits)
            cat("\nSample sizes:\n\n")
            print.data.frame(xqoi$xn, digits = digits)
          }
        cat("\n\n")
        if(!is.null(x$psclass))
          {
            qqoi <- qoi.by.sub(pscore, treat, x$psweights, x$psclass)
            cat("\nSummary of propensity score by subclasses:\n\n")
            print.data.frame(qqoi$q.table, digits = digits)
            cat("\nSample sizes by subclasses:\n\n")
            print.data.frame(qqoi$qn, digits = digits)
            cat("\n\n")
            invisible(x)
          }
      } 
  }
