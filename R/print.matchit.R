print.matchit <- function(x,  digits = max(3, getOption("digits")-3), ...)
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
    qoi <- function(xx,tt,ww,psc = NULL, full = F){
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
      if(!full){
        X.t.m <- xx[tt==1][ww[tt==1]>0]
        X.c.m <- xx[tt==0][ww[tt==0]>0]
        xsum[2,1] <- weighted.mean(X.t.m, ww[tt==1][ww[tt==1]>0])
        xsum[2,2] <- weighted.mean(X.c.m, ww[tt==0][ww[tt==0]>0])
        if(sum(tt==1)<2|(sum(tt==0)<2)){
          xsum[c(1,2),c(3,4,5)] <- NA
        } else {
          xsum[1,3] <- sd(xx,na.rm=T)
          if (is.null(x$call$exact))
            xsum[1,4] <- -1*t.test(xx~tt)$statistic
          else if(is.logical((eval(x$call$exact))))	
 		if(eval(x$call$exact))
              xsum[1,4] <- 0
            else
              xsum[1,4] <- -1*t.test(xx~tt)$statistic
          else
            xsum[1,4] <- -1*t.test(xx~tt)$statistic
          xsum[1,5] <- (mean(x1,na.rm=T)-mean(x0,na.rm=T))/sd(x1,na.rm=T)
          xsum[2,3] <- sqrt(weighted.var(xx[ww>0],ww[ww>0]))
          xsum[2,4] <- t.test.wtd(xx[ww>0],tt[ww>0],ww[ww>0])
          xsum[2,5] <- (xsum[2,1]-xsum[2,2])/sqrt(weighted.var(X.t.m,ww[tt==1][ww[tt==1]>0]))
        } 
        xn[2,1] <- sum(tt[ww>0]==1)
        xn[2,2] <- sum(tt[ww>0]==0)
        xn[2,3] <- length(tt[ww>0])
      } else {
        xn[2,] <- xn[1,] #??
        m0 <- lm(xx~ tt)
        m1 <- lm(xx~ tt + as.factor(psc))
        tdat <- m1$model
        tdat$tt <- 1
        xsum[2,1] <- mean(predict(m1,newdata=tdat)) 
        tdat$tt <- 0
        xsum[2,2] <- mean(predict(m1,newdata=tdat)) 
        xsum[1,3] <- xsum[2,3] <- sd(xx) 
        xsum[1,4] <- summary(m0)$coef[2,3]
        xsum[1,5] <- (xsum[1,1]-xsum[1,2])/sd(x1)
        xsum[2,4] <- summary(m1)$coef[2,3]
        xsum[2,5] <- (xsum[2,1]-xsum[2,2])/sd(x1)
      }
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
        if(is.null(x$match.matrix) & !identical(eval(x$call$full),TRUE))
          {
            cat("Summary of propensity score for full sample:\n\n")
            xqoi <- qoi(pscore, treat, x$psweights)
            print.data.frame(xqoi$xsum[1,], digits = digits)
            cat("\nSample sizes:\n\n")
            print.data.frame(xqoi$xn[1,], digits = digits)
          } else if(identical(eval(x$call$full),TRUE)){
            cat("\nSummary of propensity score for full and matched samples:\n\n")
            xqoi <- qoi(pscore, treat, x$psweights,
                        psc = x$data$psclass, full = T)
            print.data.frame(xqoi$xsum, digits = digits)
            cat("\nSample sizes:\n\n")
            print.data.frame(xqoi$xn, digits = digits)
          }  else {
            cat("\nSummary of propensity score for full and matched samples:\n\n")
            xqoi <- qoi(pscore, treat, x$psweights)
            print.data.frame(xqoi$xsum, digits = digits)
            cat("\nSample sizes:\n\n")
            print.data.frame(xqoi$xn, digits = digits)
          } 
        cat("\n\n")
        if(!is.null(x$psclass) & !identical(eval(x$call$full),TRUE)) {
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
