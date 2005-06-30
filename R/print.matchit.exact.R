print.matchit.exact <- function(obj){
  cat("\nCall: ", deparse(obj$call),"\n",sep = "")
  cat("\nExact Subclasses: ", max(foo$psclass),"\n",sep="")
  cat("\nSample sizes:\n")
  nn <- rbind(table(obj$treat),
              table(obj$matched,obj$treat)[2:1,])
  dimnames(nn) <- list(c("Full","Matched","Discarded"),
                       c("Control","Treated"))
  print.table(nn)
  invisible(obj)
  cat("\n")
}

