matchit.qqplot <- function(x,discrete.cutoff,
                           which.subclass=NULL, numdraws=5000,
                           interactive = T, which.xs = NULL,...){
  X <- x$X
  ## Fix X matrix so that it doesn't have any factors
  varnames <- colnames(X)
  for(var in varnames) {
        if(is.factor(X[,var])) {
                tempX <- X[,!colnames(X)%in%c(var)]
                form<-formula(substitute(~dummy-1,list(dummy=as.name(var))))
                X <- cbind(tempX, model.matrix(form, X))
        }
  }

  covariates  <- X

  if(!is.null(which.xs)){
    if(sum(which.xs%in%dimnames(covariates)[[2]])!=length(which.xs)){
      stop("which.xs is incorrectly specified")
    }
    covariates <- covariates[,which.xs,drop=F]  
  }
  treat <- x$treat
  matched <- x$weights!=0

  ratio <- x$call$ratio
  if(is.null(ratio)){ratio <- 1}
  
  ## For full or ratio matching, sample numdraws observations using the weights
  if(identical(x$call$method,"full") | (ratio!=1)) {
    t.plot <- sample(names(treat)[treat==1], numdraws/2, replace=TRUE, prob=x$weights[treat==1])
    c.plot <- sample(names(treat)[treat==0], numdraws/2, replace=TRUE, prob=x$weights[treat==0])
    
    m.covariates <- x$X[c(t.plot, c.plot),]
    m.treat <- x$treat[c(t.plot, c.plot)]
  } else {
    m.covariates <- covariates[matched,,drop=F]
    m.treat <- treat[matched]
  }
  
  if(!is.null(which.subclass)){
    subclass <- x$subclass
    sub.index <- subclass==which.subclass & !is.na(subclass)
    sub.covariates <- covariates[sub.index,,drop=F]
    sub.treat <- treat[sub.index]
    sub.matched <- matched[sub.index]
    ## Matched units in each subclass
    m.covariates <- sub.covariates[sub.matched,,drop=F]
    m.treat <- sub.treat[sub.matched]
    ## Compare to full sample--reset covariates and treat to full data set
#    covariates <- x$X
#    treat <- x$treat
  }
  nn <- dimnames(covariates)[[2]]
  nc <- length(nn)
  covariates <- data.matrix(covariates)
 # oma <- c(4, 4, 6, 4)
  oma <- c(2.25,0,3.75,1.5)
  opar <- par(mfrow = c(3, 3), mar = rep.int(1/2, 4), oma = oma)
  on.exit(par(opar))

 # par(oma=c(2.25,0,3.75,1.5))
  
  for (i in 1:nc){
    xi <- covariates[,i]
    m.xi <- m.covariates[,i]
    ni <- nn[i]
    plot(xi,type="n",axes=F)
    if(((i-1)%%3)==0){
      htext <- "QQ Plots"
      if(!is.null(which.subclass)){
        htext <- paste(htext,paste(" (Subclass ",which.subclass,")",sep=""),sep="")
      } 
      mtext(htext, 3, 2, TRUE, 0.5, cex=1.1,font=2)
      mtext("All", 3, .25, TRUE, 0.5, cex=1,font = 1)
      mtext("Matched", 3, .25, TRUE, 0.83, cex=1,font = 1)
      mtext("Control Units", 1, 0, TRUE, 2/3, cex=1,font = 1)
      mtext("Treated Units", 4, 0, TRUE, 0.5, cex=1,font = 1)
    }
    par(usr = c(0, 1, 0, 1))
    l.wid <- strwidth(nn, "user")
    cex.labels <- max(0.75, min(1.45, 0.85/max(l.wid)))
    text(0.5,0.5,ni,cex=cex.labels)
    if(length(table(xi))<=discrete.cutoff){
      xi <- jitter(xi)
      m.xi <- jitter(m.xi)
    }
    rr <- range(xi)
    eqqplot(xi[treat==0],xi[treat==1], xlim=rr,ylim=rr,axes=F,ylab="",xlab="",...)
    abline(a=0,b=1)
    abline(a=(rr[2]-rr[1])*0.1,b=1,lty=2)
    abline(a=-(rr[2]-rr[1])*0.1,b=1,lty=2)
    axis(2)
    box()
    eqqplot(m.xi[m.treat==0],m.xi[m.treat==1],xlim=rr,ylim=rr,axes=F,ylab="",xlab="",...)
    abline(a=0,b=1)
    abline(a=(rr[2]-rr[1])*0.1,b=1,lty=2)
    abline(a=-(rr[2]-rr[1])*0.1,b=1,lty=2)
    box()
    if(interactive){
      par(ask=T)
    } else {
      par(ask=F)
    }
  }
  par(ask=F)
}
