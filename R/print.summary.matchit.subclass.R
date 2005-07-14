print.summary.matchit.subclass <- function(x, digits = max(3, getOption("digits") - 3), ...){  
  sum.all <- x$sum.all
  sum.matched <- x$sum.matched
  q.table <- x$q.table
  cat("\nCall:\n", deparse(x$call),"\n\n",sep = "")
  cat("Summary of balance for all data:\n\n")
  print.data.frame(sum.all, digits=digits)
  cat("\n")
  cat("\nSummary of balance by subclasses:\n\n")
  print.table(q.table, digits = digits)
  cat("\nSample sizes by subclasses:\n\n")
  print.data.frame(x$qn, digits = digits)
  cat("\nSummary of balance across subclasses\n\n")
  print.data.frame(x$sum.subclass, digits = digits)
  cat("\nPercent Balance Improvement:\n\n")
  print.data.frame(x$reduction,digits=digits)
  cat("\n")
}
