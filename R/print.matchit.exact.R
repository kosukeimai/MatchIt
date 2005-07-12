print.matchit.exact <- function(x, digits = getOption("digits"), ...){
  cat("\nCall: ", deparse(x$call),"\n",sep = "")
  cat("\nExact Subclasses: ", max(x$subclass),"\n",sep="")
  cat("\nSample sizes:\n")
  ntab <- table(factor(x$subclass!=0,
                       levels=c("TRUE","FALSE")),
                x$treat)
  nn <- rbind(table(x$treat),
              ntab[c("TRUE","FALSE"),])
  dimnames(nn) <- list(c("Full","Matched","Discarded"),
                       c("Control","Treated"))
  print.table(nn, ...)
  invisible(x)
  cat("\n")
}

