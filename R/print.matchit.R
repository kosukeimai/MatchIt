print.matchit <- function(x, digits = getOption("digits"), ...){
  cat("\nCall: ", deparse(x$call),"\n",sep = "")
  cat("\nSample sizes:\n")
  
  #if(any(x$weights>0)) 
  #  nn <- rbind(table(x$treat),
  #              table(x$weights>0, x$treat),
  #              c(0,0))
  #else 
  #  nn <- rbind(table(x$treat),
  #              table(x$weights>0,x$treat)[2:1,])

  if(sum(x$discarded)==0) { nn <- rbind(table(x$treat), table(x$weights>0, x$treat)[2,],
		table(x$weights>0, x$treat)[1,], c(0,0)) }
  else { nn <- rbind(table(x$treat), table(x$weights>0,x$treat)[2,], 
		table(x$weights>0, x$treat)[1,], table(x$discarded, x$treat)[2,]) }
  dimnames(nn) <- list(c("Full","Matched","Unmatched","Discarded"),
                       c("Control","Treated"))
  print.table(nn, ...)
  invisible(x)
  cat("\n")
}
