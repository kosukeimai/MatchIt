print.matchit <- function(obj){
  cat("\nCall: ", deparse(obj$call),"\n",sep = "")
  cat("\nSample sizes:\n")
  nn <- rbind(table(obj$data$treat),
              table(obj$matched,obj$data$treat)[2:1,])
  dimnames(nn) <- list(c("Full","Matched","Discarded"),
                       c("Control","Treated"))
  print.table(nn)
  invisible(obj)
  cat("\n")
}
