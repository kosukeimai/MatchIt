hist.pscore <- function(x, numdraws=5000){
  treat <- x$treat
  pscore <- x$distance
  weights <- x$weights
  matched <- weights!=0
  q.cut <- x$q.cut
  cwt <- sqrt(weights)
  ratio <- x$call$ratio
  if(is.null(ratio)){ratio <- 1}
  
  ## For full or ratio matching, sample numdraws observations using the weights
  if(identical(x$call$method,"full") | (ratio!=1)) {
    pscore.treated.matched <- sample(names(treat)[treat==1],
                                     numdraws/2, replace=TRUE,
                                     prob=x$weights[treat==1])
    pscore.treated.matched <- pscore[pscore.treated.matched]
    pscore.control.matched <- sample(names(treat)[treat==0],
                                     numdraws/2, replace=TRUE,
                                     prob=x$weights[treat==0])
    pscore.control.matched <- pscore[pscore.control.matched]
    
  } else {
    pscore.treated.matched <- pscore[treat==1 & weights!=0]
    pscore.control.matched <- pscore[treat==0 & weights!=0]
  }
  par(mfrow=c(2,2))
  hist(pscore[treat==1],xlim=range(na.omit(pscore)),
       xlab="Propensity Score", freq=F,
       main="Raw Treated")
  hist(pscore.treated.matched,xlim=range(na.omit(pscore)),
       xlab="Propensity Score", freq=F,
       main="Matched Treated")
  if(!is.null(q.cut)){abline(v=q.cut,col="grey",lty=1)}
  hist(pscore[treat==0],xlim=range(na.omit(pscore)),
       xlab="Propensity Score", freq=F,
       main="Raw Control")
  hist(pscore.control.matched,xlim=range(na.omit(pscore)),
       xlab="Propensity Score", freq=F,
       main="Matched Control")
  if(!is.null(q.cut)){abline(v=q.cut,col="grey",lty=1)}
}
