summary.matchit <- function(object, verbose=F, sig=2, ...) {
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
  qoi <- function(xx,tt,ww, psc = NULL, full = T){
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
        xsum[1,4] <- -1*t.test(xx~tt)$sta
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
      xsum[1,3] <- xsum[2,3] <- sd(foo1$data$pscore) 
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
  #setting up covariates
  pscore <- object$data[,"pscore"]
  treat <- eval(object$treat,object$data)
  names(treat) <- names(pscore)
  data <- object$data
  weights <- object$psweights
  covariates <- as.data.frame(model.matrix(terms(object$formula),data)[,,drop=FALSE])
  if ("(Intercept)"%in%colnames(covariates))
    covariates <- covariates[,-match("(Intercept)", colnames(covariates)),drop=FALSE]
  mahvars <- eval(object$call$mahvars)
  exact <- eval(object$call$exact)
  subclass <- eval(object$call$subclass)
  if(is.null(subclass)){subclass <- 0}
  psclass <- object$psclass
  if (!is.null(mahvars)){
    w <- mahvars%in%names(covariates)
    if(sum(w)!=length(mahvars)){
      md <- as.data.frame(as.matrix(data[,mahvars[!w]]))
      names(md) <- mahvars[!w]
      row.names(md) <- row.names(covariates)
      covariates <- cbind(covariates,md)
    }
  }
  if (!identical(TRUE, exact) & !identical(FALSE, exact) & !is.null(exact)) {
    w <- eval(exact)%in%names(covariates)
    if(sum(w)!=length(exact)) {
      ed <- as.data.frame(as.matrix(data[,exact[!w]]))
      names(ed) <- exact[!w]
      row.names(ed) <- row.names(covariates)
      covariates <- cbind(covariates,ed)
    }    
  }
  if(!identical(TRUE,exact)){ 
    covariates<- cbind(pscore,covariates)
  } 

  if(!identical(eval(object$call$full),TRUE)){
    aa <- apply(covariates,2,qoi,tt=treat,ww=weights)
  } else {
    aa <- apply(covariates,2,qoi,tt=treat,ww=weights,
                psc = object$data$psclass, full = T)
  }
  k <- length(aa)
  sum.all <- as.data.frame(matrix(0,k,5))
  sum.matched <- as.data.frame(matrix(0,k,6))
  row.names(sum.all) <- row.names(sum.matched) <- names(aa)
  names(sum.all) <- names(sum.matched) <- names(aa[[1]]$xsum)
  names(sum.matched) <- c(names(aa[[1]]$xsum),"Reduction")
  sum.all.int <- sum.matched.int <- NULL
  for(i in 1:k){
    sum.all[i,] <- aa[[i]]$xsum[1,]
    sum.matched[i,1:5] <- aa[[i]]$xsum[2,]
    sum.matched[i,6] <- (abs(sum.all[i,5])-abs(sum.matched[i,5]))>0
    if(verbose){
      for(j in i:k){
        x2 <- covariates[,i]*as.matrix(covariates[,j])
        sum.all.int <- rbind(sum.all.int,qoi(x2,tt=treat,ww=weights)$xsum[1,])
        Reduction <- (abs(sum.all.int[nrow(sum.all.int),5])-abs(qoi(x2,tt=treat,ww=weights)$xsum[2,5]))>0
        sum.matched.int <- rbind(sum.matched.int,cbind(qoi(x2,tt=treat,ww=weights)$xsum[2,],Reduction))
        row.names(sum.all.int)[nrow(sum.all.int)] <-
          row.names(sum.matched.int)[nrow(sum.matched.int)] <-
            paste(names(covariates)[i],names(covariates)[j],sep="x")
      }
    }
  }
  xn <- aa[[1]]$xn
  sum.all <- rbind(sum.all,sum.all.int)
  sum.matched <- rbind(sum.matched,sum.matched.int)
  
  #now subclassification
  if(!identical(subclass,0) & !identical(eval(object$call$full),TRUE))
    {
      qbins <- max(psclass,na.rm=TRUE)
      if(verbose){
        q.table <- array(0,dim=c(k+sum(1:k),5,qbins))
        ii <- 0
        nn <- NULL
      } else {
        q.table <- array(0,dim=c(k,5,qbins))
      }
      aa <- apply(covariates,2,qoi.by.sub,tt=treat,ww=weights,
                  qq=object$psclass)
      for(i in 1:k){
        if(!verbose){
          q.table[i,,] <- as.matrix(aa[[i]]$q.table)
          nn <- names(aa)
        } else {
          ii <- ii + 1 
          q.table[ii,,] <- as.matrix(aa[[i]]$q.table)
          nn <- c(nn,names(aa)[i])
          for(j in i:k){
            ii <- ii + 1 
            x2 <- covariates[,i]*as.matrix(covariates[,j])
            q.table[ii,,] <- as.matrix(qoi.by.sub(x2,tt=treat,ww=weights,qq=object$psclass)$q.table)
            nn <- c(nn,paste(names(covariates)[i],names(covariates)[j],sep="x"))
          }
        }
      }   
      qn <- aa[[1]]$qn
      dimnames(q.table) <- list(nn,row.names(aa[[i]]$q.table),paste("Subclass",1:qbins))
    } else {
      q.table <- NULL; q.bias <- NULL; qn <- NULL
    }
  if(identical(exact,TRUE)){
    qbins <- max(object$psclass,na.rm=TRUE)
    q.table <- as.data.frame(matrix(0,qbins,k+3))
    names(q.table) <- c(names(covariates),"Treated","Control","Total")
    for(i in 1:qbins){
      qi <- object$psclass==i
      q.table[i,] <- c(as.numeric(covariates[qi,,drop=F][1,]), sum(treat[qi]==1), sum(treat[qi]==0), length(treat[qi]))
    }
  }
  ans <- list(sum.all = sum.all, sum.matched = sum.matched,
              q.table=q.table, call=object$call, verbose=verbose,
              xn = xn, match.matrix = object$match.matrix, sig = sig,
              psclass=object$psclass, in.sample = object$in.sample,
              qn = qn)
  class(ans) <- "summary.matchit"
  ans
}  

