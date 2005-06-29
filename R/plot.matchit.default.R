plot.matchit.default <- function(obj,discrete.cutoff=5){
  match.matrix <- obj$match.matrix
  xdata <- obj$data
  treata <- model.frame(obj$formula,xdata)[,1,drop=FALSE]
  pscore <- xdata[,"pscore"]
  in.sample <- obj$in.sample
  covariates <- obj$covariates
  treat <- as.vector(treata[,1])
  names(treat) <- row.names(treata)
  covariates <- model.frame(delete.response(terms(obj$formula)),xdata)[,,drop=FALSE]
  mahvars <- eval(obj$call$mahvars)
  exact <- eval(obj$call$exact)
  if (!is.null(mahvars)){
    w <- mahvars%in%names(covariates)
    if(sum(w)!=length(mahvars)){
      md <- as.data.frame(as.matrix(xdata[,mahvars[!w]]))
      names(md) <- mahvars[!w]
      covariates <- cbind(covariates,md)
    }
  }
  if (!identical(TRUE, exact) & !identical(FALSE, exact)) {
    w <- exact%in%names(covariates)
    if(sum(w)!=length(exact)) {
      ed <- as.data.frame(as.matrix(xdata[,exact[!w]]))
      names(ed) <- exact[!w]
      covariates <- cbind(covariates,ed)
    }
  }
  nn <- names(covariates)
  nc <- length(nn)
  covariates <- data.matrix(covariates)
  oma <- c(4, 4, 6, 4)
  opar <- par(mfrow = c(3, 3), mar = rep.int(1/2, 4), oma = oma, ask=T)
  on.exit(par(opar))
  for (i in 1:nc){
    xi <- covariates[,i]
    ni <- nn[i]
    plot(xi,type="n",axes=F)
    if(((i-1)%%3)==0){
      mtext("QQ Plots", 3, 3, TRUE, 0.5, cex=1.2,font=2)
      mtext("Raw", 3, 1, TRUE, 0.5, cex=1.2,font = 1)
      mtext("Matched", 3, 1, TRUE, 0.83, cex=1.2,font = 1)
      mtext("Control Units", 1, 0, TRUE, 2/3, cex=1,font = 1)
      mtext("Treated Units", 4, 0, TRUE, 0.5, cex=1,font = 1)
    }
    par(usr = c(0, 1, 0, 1))
    l.wid <- strwidth(nn, "user")
    cex.labels <- max(0.8, min(2, 0.9/max(l.wid)))
    text(0.5,0.5,ni,cex=cex.labels)
    if(length(table(xi))>discrete.cutoff){
      rr <- range(xi)
      qqplot(xi[treat==0],xi[treat==1],
             xlim=rr,ylim=rr,axes=F,ylab="",
             xlab="")
      abline(a=0,b=1)
      abline(a=(rr[2]-rr[1])*0.1,b=1,lty=2)
      abline(a=-(rr[2]-rr[1])*0.1,b=1,lty=2)
      axis(2)
      box()
      qqplot(xi[treat==0 & obj$matched],
             xi[treat==1 & obj$matched],
             xlim=rr,ylim=rr,axes=F,ylab="",
             xlab="")
      abline(a=0,b=1)
      abline(a=(rr[2]-rr[1])*0.1,b=1,lty=2)
      abline(a=-(rr[2]-rr[1])*0.1,b=1,lty=2)
      box()
    } else{
      tb1 <- table(xi,treat)
      tb1 <- t(t(tb1)/apply(tb1,2,sum))
      tb2 <- table(xi,treat,obj$matched)[,,"TRUE"]
      tb2 <- t(t(tb2)/apply(tb2,2,sum))
      rr <- range(tb1,tb2)
      qqplot(xi[treat==0],xi[treat==1],
             xlim=rr,ylim=rr,axes=F,
             type="n",xlab="",ylab="")
      abline(a=0,b=1)
      abline(a=(rr[2]-rr[1])*0.1,b=1,lty=2)
      abline(a=-(rr[2]-rr[1])*0.1,b=1,lty=2)
      axis(2)
      text(tb1[,1],tb1[,2],row.names(tb1))
      box()
      qqplot(xi[treat==0 & obj$matched],
             xi[treat==1 & obj$matched],
             xlim=rr,ylim=rr,axes=F,
             type="n",xlab="",ylab="")
      text(tb2[,1],tb2[,2],row.names(tb2))
      abline(a=0,b=1)
      abline(a=(rr[2]-rr[1])*0.1,b=1,lty=2)
      abline(a=-(rr[2]-rr[1])*0.1,b=1,lty=2)
      box()
    }
  }  
}
