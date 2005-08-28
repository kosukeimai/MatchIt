print.summary.matchit.exact <- function(x, digits = max(3,
                                             getOption("digits") - 3),
                                        ...){  
  cat("\nCall:", deparse(x$call), sep = "\n")
  cat("\nSample sizes:\n")
  ntab <- table(factor(!is.na(x$subclass),
                      levels=c("TRUE","FALSE")), x$treat)
  nn <- rbind(table(x$treat),
             ntab[c("TRUE","FALSE"),])
  dimnames(nn) <- list(c("All","Matched","Discarded"),
                       c("Control","Treated"))
  print.table(nn,digits=digits)
  cat("\nMatched sample sizes by subclass:\n")
  print.data.frame(x$q.table, digits = digits)
  cat("\n")
  invisible(x)
}
