print.summary.matchit.exact <- function(obj, digits = max(3, getOption("digits") - 3), ...){  
  cat("\nSample sizes for full and exactly matched data:\n\n")
  ntab <- table(factor(obj$subclass!=0,
                       levels=c("TRUE","FALSE")),
                obj$treat)
  nn <- rbind(table(obj$treat),
              ntab[c("TRUE","FALSE"),])
  dimnames(nn) <- list(c("Full","Matched","Discarded"),
                       c("Control","Treated"))
  print.table(nn,digits=digits)
  cat("\n")
  if(!obj$verbose){
    cat("\nMatched sample sizes by subclass:\n\n")
  } else {
    cat("\nMatched sample sizes by subclass with covariates:\n\n")
  }
  print.data.frame(obj$q.table, digits = digits)
  cat("\n")
}
