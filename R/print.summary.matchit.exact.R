print.summary.matchit.exact <- function(x, digits = max(3, getOption("digits") - 3), ...){  
  #cat("\nSample sizes for full and exactly matched data:\n\n")
  #ntab <- table(factor(!is.na(x$subclass),
  #                     levels=c("TRUE","FALSE")), x$treat)
  #nn <- rbind(table(x$treat),
  #            ntab[c("TRUE","FALSE"),])
  #dimnames(nn) <- list(c("Full","Matched","Discarded"),
  #                     c("Control","Treated"))
  #print.table(nn,digits=digits)

  cat("\nCall:\n", deparse(x$call),"\n\n",sep = "")
  cat("\nMatched sample sizes by subclass:\n\n")
  print.data.frame(x$q.table, digits = digits)
  cat("\n")
  invisible(x)
}
