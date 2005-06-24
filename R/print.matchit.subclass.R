print.matchit.subclass <- function(obj){
  cat("\nCall: ", deparse(obj$call),"\n",sep = "")
  cat("\nSample sizes by subclasses:\n\n")
  nsub <- table(obj$psclass,obj$data$treat)
  nn <- rbind(table(obj$data$treat),nsub)
  dimnames(nn) <-
    list(c("Full",paste("Subclass",dimnames(nsub)[[1]])),
         c("Control","Treated"))
  print.table(nn)
  invisible(obj)
  cat("\n")
}
