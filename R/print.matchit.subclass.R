print.matchit.subclass <- function(x, digits = getOption("digits"), ...){
  cat("\nCall: ", deparse(x$call), sep = "\n")
  cat("\nSample sizes by subclasses:\n\n")
  nsub <- table(x$subclass,x$treat)
  nn <- rbind(table(x$treat),nsub)
  dimnames(nn) <-
    list(c("All",paste("Subclass",dimnames(nsub)[[1]])),
         c("Control","Treated"))
  print.table(nn, ...)
  invisible(x)
  cat("\n")
}
