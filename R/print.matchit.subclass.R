print.matchit.subclass <- function(obj){
  cat("\nCall: ", deparse(obj$call),"\n",sep = "")
  cat("\nSample sizes by subclasses:\n\n")
  nsub <- table(obj$subclass,obj$treat)
  nn <- rbind(table(obj$treat),nsub)
  dimnames(nn) <-
    list(c("Full",paste("Subclass",dimnames(nsub)[[1]])),
         c("Control","Treated"))
  print.table(nn)
  invisible(obj)
  cat("\n")
}
