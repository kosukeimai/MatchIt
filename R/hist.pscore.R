hist.pscore <- function(x){
  treat <- x$treat
  pscore <- x$distance
  weights <- x$weights
  matched <- weights!=0
  q.cut <- x$q.cut
  cwt <- sqrt(weights)
  par(mfrow=c(2,2))
  hist(pscore[treat==1],xlim=range(na.omit(pscore)),
       xlab="Propensity Score", freq=F,
       main="Raw Treated")
  hist(pscore[treat==0],xlim=range(na.omit(pscore)),
       xlab="Propensity Score", freq=F,
       main="Raw Control")
  hist(pscore[treat==1 & weights!=0],xlim=range(na.omit(pscore)),
       xlab="Propensity Score", freq=F,
       main="Matched Treated")
  hist(pscore[treat==0 & weights!=0],xlim=range(na.omit(pscore)),
       xlab="Propensity Score", freq=F,
       main="Matched Control")  
}
