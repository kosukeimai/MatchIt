print.matchit <- function(obj){
  cat("\nCall: ", deparse(obj$call),"\n",sep = "")
  cat("\nSample sizes:\n")
  nn <- rbind(table(obj$treat),
              table(obj$weights!=0,obj$treat)[2:1,])
  dimnames(nn) <- list(c("Full","Matched","Discarded"),
                       c("Control","Treated"))
  print.table(nn)
  invisible(obj)
  cat("\n")
}
