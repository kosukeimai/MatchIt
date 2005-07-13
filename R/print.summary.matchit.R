print.summary.matchit <- function(x, digits = max(3, getOption("digits") - 3), ...){  

  sum.all <- x$sum.all
  sum.matched <- x$sum.matched
  q.table <- x$q.table
  xn <- x$nn
  cat("\nCall:\n", deparse(x$call),"\n\n",sep = "")
  cat("Summary of balance for all data:\n\n")
  print.data.frame(sum.all, digits=digits)
  cat("\n")

  xs1 <- sum.matched
  cc <- row.names(sum.all)
  if(!is.null(x$match.matrix) | identical(eval(x$call$full),TRUE))
    {
      cat("Summary of balance for matched data:\n\n")
      print.data.frame(xs1, digits=digits)
      cat("\nPercent Balance Improvement:\n\n")
      print.data.frame(x$reduction,digits=digits)
      cat("\nSample sizes:\n\n")
      print.table(xn, digits=digits)
      cat("\n")
    }
  invisible(x)
}
