print.summary.matchit <- function(x, digits = max(3, getOption("digits") - 3), ...){  
  verbose <- x$verbose
  sum.all <- x$sum.all
  sum.matched <- x$sum.matched
  q.table <- x$q.table
  xn <- x$xn
  sig <- x$sig
  cat("\nAssignment model specification:\n", deparse(x$call),"\n\n",sep = "")
  if(verbose){
    cat("Summary of covariates and interactions for all data:\n\n")
  } else{
    cat("Summary of covariates for all data:\n\n")
  }
  print(sum.all, digits=digits)
  cat("\n")
  xs1 <- sum.matched
  cc <- row.names(sum.all)
  if(!is.null(x$match.matrix))
    {
      if(verbose){
        cat("Summary of covariates and interactions for matched data:\n\n")
      } else {
        cat("Summary of covariates for matched data:\n\n")
      }
      print(xs1, digits=digits)
      cat("\nSample sizes:\n\n")
      print.data.frame(xn, digits=digits)
      cat("\n")
      cat("Problematic covariates:  ")
      xs2 <- xs1[1:(nrow(xs1)-1),]
      cat(row.names(xs2)[!is.na(xs2[,4]) & abs(xs2[,4])>sig])
      cat("\n")
    }
  if(!is.null(x$psclass) & !identical(eval(x$call$opt),TRUE))
    {
      if(identical(eval(x$call$exact),TRUE)){
        cat("\nSample sizes for full and exactly matched data:\n\n")
        print(xn,digits = digits)
        cat("\nSample sizes by covariates:\n\n")
        print(q.table, digits = digits)
        cat("\n")
      } else {
        print(q.table, digits = digits)
        cat("\nSample sizes by subclasses:\n\n")
        print.data.frame(x$qn, digits = digits)
        cat("\n")
        cat("Problematic covariates:\n")
        for(i in 1:dim(q.table)[3])
          {
            cat("Subclass ", i, ":  ")
            xs1 <- q.table[cc,,i]
            cat(row.names(xs1)[!is.na(xs1[,4]) & abs(xs1[,4])>sig])
            cat("\n")
          }
      }
    }
  
  cat("Number of units discarded:  ", sum(x$in.sample==FALSE))
  cat("\n\n")
  invisible(x)
}
