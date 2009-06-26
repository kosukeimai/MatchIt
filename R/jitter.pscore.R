jitter.pscore <- function(x, interactive,pch=1,cex=NULL,...){
  treat <- x$treat
  pscore <- x$distance
  weights <- x$weights
  matched <- weights!=0
  q.cut <- x$q.cut
  jitp <- jitter(rep(1,length(treat)),factor=6)+(treat==1)*(weights==0)-(treat==0)-(weights==0)*(treat==0)
  cwt <- sqrt(weights)
  minp <- min(pscore,na.rm=T)
  maxp <- max(pscore,na.rm=T)
  plot(pscore,xlim=c(minp,maxp+0.1*(maxp-minp)),ylim=c(-1.5,2.5),
       type="n",ylab="",xlab="Propensity Score",
       axes=F,main="Distribution of Propensity Scores",...)
  if(!is.null(q.cut)){abline(v=q.cut,col="grey",lty=1)}
  if(is.null(cex)){
   points(pscore[treat==1&weights!=0],jitp[treat==1&weights!=0],
        pch=pch,cex=cwt[treat==1&weights!=0],...)
   points(pscore[treat==0&weights!=0],jitp[treat==0&weights!=0],
        pch=pch,cex=cwt[treat==0&weights!=0],...)
   points(pscore[treat==1&weights==0],jitp[treat==1&weights==0],
        pch=pch,cex=1,...)
   points(pscore[treat==0&weights==0],jitp[treat==0&weights==0],
        pch=pch,cex=1,...)
  }else{
    points(pscore[treat==1&weights!=0],jitp[treat==1&weights!=0],
        pch=pch,cex=cex,...)
    points(pscore[treat==0&weights!=0],jitp[treat==0&weights!=0],
         pch=pch,cex=cex,...)
    points(pscore[treat==1&weights==0],jitp[treat==1&weights==0],
         pch=pch,cex=cex,...)
    points(pscore[treat==0&weights==0],jitp[treat==0&weights==0],
         pch=pch,cex=cex,...)
  }
  axis(1)
  text(sum(range(na.omit(pscore)))/2,2.5,"Unmatched Treatment Units")
  text(sum(range(na.omit(pscore)))/2,1.5,"Matched Treatment Units")
  text(sum(range(na.omit(pscore)))/2,0.5,"Matched Control Units")
  text(sum(range(na.omit(pscore)))/2,-0.5,"Unmatched Control Units")
  box()
  if(interactive==TRUE) {
	print("To identify the units, use first mouse button; to stop, use second.")
  	identify(pscore,jitp,names(treat),atpen=T)
  }
}
