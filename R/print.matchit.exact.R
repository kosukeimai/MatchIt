print.matchit.exact <- function(obj){
  cat("\nCall: ", deparse(obj$call),"\n",sep = "")
  cat("\nExact Subclasses: ", max(obj$subclass),"\n",sep="")
  cat("\nSample sizes:\n")
  ntab <- table(factor(obj$subclass!=0,
                       levels=c("TRUE","FALSE")),
                obj$treat)
  nn <- rbind(table(obj$treat),
              ntab[c("TRUE","FALSE"),])
  dimnames(nn) <- list(c("Full","Matched","Discarded"),
                       c("Control","Treated"))
  print.table(nn)
  invisible(obj)
  cat("\n")
}

