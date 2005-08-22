print.matchit <- function(x, digits = getOption("digits"), ...){
  cat("\nCall: ", deparse(x$call),"\n", sep="\n")
  cat("\nSample sizes:\n")
  
  #if(any(x$weights>0)) 
  #  nn <- rbind(table(x$treat),
  #              table(x$weights>0, x$treat),
  #              c(0,0))
  #else 
  #  nn <- rbind(table(x$treat),
  #              table(x$weights>0,x$treat)[2:1,])

  nn <- matrix(0, ncol=2, nrow=4)
  nn[1,] <- c(sum(x$treat==0), sum(x$treat==1))
  nn[2,] <- c(sum(x$treat==0 & x$weights>0), sum(x$treat==1 & x$weights>0))
  nn[3,] <- c(sum(x$treat==0 & x$weights==0 & x$discarded==0), sum(x$treat==1 & x$weights==0 & x$discarded==0))
  nn[4,] <- c(sum(x$treat==0 & x$weights==0 & x$discarded==1), sum(x$treat==1 & x$weights==0 & x$discarded==1))

  dimnames(nn) <- list(c("Full","Matched","Unmatched","Discarded"),
                       c("Control","Treated"))
  print.table(nn, ...)
  invisible(x)
  cat("\n")
}
