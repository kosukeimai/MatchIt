summary.neyman<-function(object,...){
  tval<-object$ate/object$se
  res<-cbind(object$ate, object$se, tval, 2*pnorm(abs(tval), lower.tail=FALSE))
  colnames(res)<-c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  rownames(res)[1]<-c("Overall")
  
  sampsize <- cbind(object$Ntrt, object$Ncont)
  colnames(sampsize) <- c("Ntrt", "Ncont")
  rownames(sampsize) <- rownames(res)

  output <- list(results=res, sample.size=sampsize)
  class(output)<-"summary.neyman"
  output
}
