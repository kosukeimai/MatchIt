print.neyman<-function(x, digits=max(3, getOption("digits") -3), ...){
  print.res<-cbind(x$ate, x$se)
  colnames(print.res)<-c("ATE", "SD")
  print.res <- print.res[1,]
                                                                                                                                                             
  cat("\n Average Treatment Effect:\n \n")
  print.matrix(print.res, digits=digits, ...)
  cat("\n")
}
