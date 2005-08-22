print.summary.matchit <- function(x, digits = max(3,
                                       getOption("digits") - 3), ...){
  
  sum.all <- x$sum.all
  sum.matched <- x$sum.matched
  q.table <- x$q.table
  xn <- x$nn
  cat("\nCall:\n", deparse(x$call),"\n",sep = "\n")
  cat("Summary of balance for all data:\n\n")
  print.data.frame(round(sum.all,digits))
  cat("\n")

  xs1 <- sum.matched
  cc <- row.names(sum.all)
  if(!is.null(x$sum.matched) | identical(eval(x$call$method),"full"))
    {
      cat("Summary of balance for matched data:\n\n")
      print.data.frame(round(xs1,digits))
      cat("\nPercent Balance Improvement:\n\n")
      print.data.frame(round(x$reduction,digits))
      cat("\nSample sizes:\n\n")
      print.table(xn, digits=digits)
      cat("\n")
    }
  invisible(x)
}
