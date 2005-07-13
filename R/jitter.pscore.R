jitter.pscore <- function(obj){
  treat <- obj$treat
  pscore <- obj$distance
  weights <- obj$weights
  matched <- weights!=0
  q.cut <- obj$q.cut
  jitp <- jitter(rep(1,length(treat)),factor=6)-(treat==0)
  cwt <- sqrt(weights)
  plot(pscore,xlim=range(na.omit(pscore)),ylim=c(-1,2),
       type="n",ylab="",xlab="Propensity Score",
       axes=F,main="Distribution of Propensity Scores")
  if(!is.null(q.cut)){abline(v=q.cut,col="grey",lty=1)}
  points(pscore[treat==1&weights!=0],jitp[treat==1&weights!=0],
         pch=18,cex=cwt[treat==1&weights!=0])
  points(pscore[treat==1&weights==0],jitp[treat==1&weights==0],
         pch=5,col="grey",cex=0.5)
  points(pscore[treat==0&weights!=0],jitp[treat==0&weights!=0],
         pch=18,cex=cwt[treat==0&weights!=0])
  points(pscore[treat==0&weights==0],jitp[treat==0&weights==0],
         pch=5,col="grey",cex=0.5)
  axis(1)
  text(sum(range(na.omit(pscore)))/2,1.5,"Treatment Units")
  text(sum(range(na.omit(pscore)))/2,-0.5,"Control Units")
  box()
  print("To identify the units, use first mouse button; to stop, use second.")
  identify(pscore,jitp,names(treat))
}
