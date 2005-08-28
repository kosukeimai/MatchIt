print.matchit.exact <- function(x, digits = getOption("digits"), ...){
  cat("\nCall: ", deparse(x$call), sep = "\n")
  cat("\nExact Subclasses: ", max(x$subclass, na.rm=T),"\n",sep="")
  cat("\nSample sizes:\n")
  ntab <- table(factor(!is.na(x$subclass),
                       levels=c("TRUE","FALSE")),
                x$treat)
  nn <- rbind(table(x$treat),
              ntab[c("TRUE","FALSE"),])
  dimnames(nn) <- list(c("All","Matched","Unmatched"),
                       c("Control","Treated"))
  print.table(nn, ...)
  invisible(x)
  cat("\n")
}

