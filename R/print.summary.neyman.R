print.summary.neyman<-function(x, digits=max(3, getOption("digits")-3),
                               signif.stars=getOption("show.signif.stars"), 
                               ...){
  cat("\n Average Treatment Effect:\n")
  printCoefmat(x$results, digits=digits, signif.starts = signif.stars, ...)
  cat("\n")

  cat("\n Sample Sizes:\n")
  print.matrix(x$sample.size)
  cat("\n")
}
