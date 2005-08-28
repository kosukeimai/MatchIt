print.summary.matchit.subclass <- function(x, digits = max(3,
                                                getOption("digits") -
                                                3), ...){   
  sum.all <- x$sum.all
  sum.matched <- x$sum.matched
  q.table <- x$q.table
  cat("\nCall:", deparse(x$call), sep = "\n")
  cat("Summary of balance for all data:\n")
  print.data.frame(round(sum.all,digits))
  cat("\n")
  cat("\nSummary of balance by subclasses:\n")
  print.table(round(q.table, digits))
  cat("\nSample sizes by subclasses:\n")
  print.data.frame(x$qn, digits = digits)
  cat("\nSummary of balance across subclasses\n")
  print.data.frame(round(x$sum.subclass, digits))
  cat("\nPercent Balance Improvement:\n")
  print.data.frame(round(x$reduction,digits))
  cat("\n")
}
