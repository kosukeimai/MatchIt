print.summary.matchit <- function(x, digits = max(3, getOption("digits") - 3), ...){  
  verbose <- x$verbose
  sum.all <- x$sum.all
  sum.matched <- x$sum.matched
  q.table <- x$q.table
  xn <- x$nn
  cat("\nAssignment model specification:\n", deparse(x$call),"\n\n",sep = "")
  if(verbose){
    cat("Summary of covariates and interactions for all data:\n\n")
  } else{
    cat("Summary of covariates for all data:\n\n")
  }
  print.data.frame(sum.all, digits=digits)
  cat("\n")
  xs1 <- sum.matched
  cc <- row.names(sum.all)
  if(!is.null(x$match.matrix) | identical(eval(x$call$full),TRUE))
    {
      if(verbose){
        cat("Summary of covariates and interactions for matched data:\n\n")
      } else {
        cat("Summary of covariates for matched data:\n\n")
      }
      print.data.frame(xs1, digits=digits)
      cat("\nPercent Balance Improvement:\n\n")
      print.data.frame(x$reduction,digits=digits)
      cat("\nSample sizes:\n\n")
      print.table(xn, digits=digits)
      cat("\n")
    }
  if(!is.null(x$psclass) & !identical(eval(x$call$full),TRUE))
    {
      if(identical(eval(x$call$exact),TRUE)){
        cat("\nSample sizes for full and exactly matched data:\n\n")
        print.table(xn,digits = digits)
        cat("\nSample sizes by covariates:\n\n")
        print.table(q.table, digits = digits)
        cat("\n")
      } else {
        print.table(q.table, digits = digits)
        cat("\nSample sizes by subclasses:\n\n")
        print.data.frame(x$qn, digits = digits)
        cat("\n")
      }
    }
  invisible(x)
}
