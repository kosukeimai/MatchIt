print.summary.matchit.subclass <- function(x, digits = max(3, getOption("digits") - 3), ...){  
  verbose <- x$verbose
  sum.all <- x$sum.all
  sum.matched <- x$sum.matched
  q.table <- x$q.table
  cat("\nAssignment model specification:\n", deparse(x$call),"\n\n",sep = "")
  if(verbose){
    cat("Summary of covariates and interactions for all data:\n\n")
  } else{
    cat("Summary of covariates for all data:\n\n")
  }
  print.data.frame(sum.all, digits=digits)
  cat("\n")
  cat("\nSummary of covariates by subclasses:\n\n")
  print.table(q.table, digits = digits)
  cat("\nSample sizes by subclasses:\n\n")
  print.data.frame(x$qn, digits = digits)
  cat("\n")
}
