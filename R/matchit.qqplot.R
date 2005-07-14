matchit.qqplot <- function(x,discrete.cutoff,which.subclass=NULL){
  covariates  <- x$X
  treat <- x$treat
  matched <- x$weights!=0
  if(!is.null(which.subclass)){
    subclass <- x$subclass
    sub.index <- subclass==which.subclass & !is.na(subclass)
    covariates <- covariates[sub.index,]
    treat <- treat[sub.index]
    matched <- matched[sub.index]
  }
  nn <- dimnames(covariates)[[2]]
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
      htext <- "QQ Plots"
      if(!is.null(which.subclass)){
        htext <- paste(htext,paste(" (Subclass ",which.subclass,")",sep=""),sep="")
      } 
      mtext(htext, 3, 3, TRUE, 0.5, cex=1.2,font=2)
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
      qqplot(xi[treat==0 & matched],
             xi[treat==1 & matched],
             xlim=rr,ylim=rr,axes=F,ylab="",
             xlab="")
      abline(a=0,b=1)
      abline(a=(rr[2]-rr[1])*0.1,b=1,lty=2)
      abline(a=-(rr[2]-rr[1])*0.1,b=1,lty=2)
      box()
    } else{
      xij <- jitter(xi)
      rr <- range(xij)
      qqplot(xij[treat==0],xij[treat==1],
             xlim=rr,ylim=rr,axes=F,
             xlab="",ylab="")
      abline(a=0,b=1)
      abline(a=(rr[2]-rr[1])*0.1,b=1,lty=2)
      abline(a=-(rr[2]-rr[1])*0.1,b=1,lty=2)
      axis(2)
      box()
      qqplot(xij[treat==0 & matched],
             xij[treat==1 & matched],
             xlim=rr,ylim=rr,axes=F,
             xlab="",ylab="")
      abline(a=0,b=1)
      abline(a=(rr[2]-rr[1])*0.1,b=1,lty=2)
      abline(a=-(rr[2]-rr[1])*0.1,b=1,lty=2)
      box()
    }
  }  
}
