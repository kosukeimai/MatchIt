print.summary.matchit <- function(x, digits = max(3,
                                       getOption("digits") - 3), ...){
  
  sum.all <- x$sum.all
  sum.matched <- x$sum.matched
  q.table <- x$q.table
  xn <- x$nn
  cat("\nCall:", deparse(x$call), sep = "\n")
  cat("\nSummary of balance for all data:\n")
  print.data.frame(round(sum.all,digits))
  cat("\n")

  xs1 <- sum.matched
  cc <- row.names(sum.all)
  if(!is.null(x$sum.matched) | identical(eval(x$call$method),"All"))
    {
      cat("\nSummary of balance for matched data:\n")
      print.data.frame(round(xs1,digits))
      cat("\nPercent Balance Improvement:\n")
      print.data.frame(round(x$reduction,digits))
      cat("\nSample sizes:\n")
      print.table(xn, digits=digits)
      cat("\n")
    }
  invisible(x)
}
