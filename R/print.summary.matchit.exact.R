print.summary.matchit.exact <- function(x, digits = max(3,
                                             getOption("digits") - 3),
                                        ...){  
  cat("\nCall:", deparse(x$call), sep = "\n")
  cat("\nSample sizes:\n")
  print.table(x$nn,digits=digits)
  cat("\nMatched sample sizes by subclass:\n")
  print.data.frame(x$q.table, digits = digits)
  cat("\n")
  invisible(x)
}
