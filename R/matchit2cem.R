# matchit2cem - matchit wrapper for cem matching algorithm
#
# 06/10/2008 - m.blackwell
#
# this function takes inputs from matchit() and returns the
# strata for each observation in the subclass entry and the
# weight for each observation in the weight entry. No match
# matrix is returned since matches are not unique within
# strata. 
#
matchit2cem <- function(treat, X, data, distance, discarded, is.full.mahalanobis,
                            ratio = 1, verbose = FALSE, k2k.method=NULL, ...) {

  
  #if (!("cem" %in% .packages(all = TRUE)))
  #  install.packages("cem",repos="http://gking.harvard.edu/")
  requireNamespace(cem)

  if (verbose)
    cat("Coarsened exact matching...\n")
  
  n <- length(treat)
  
  # cem takes the data all together and wants the treatment specified
  # with the column name of the data frame. Here we massage the matchit
  # inputs to this format. Note that X has its proper columnames, but
  # treat does not have the original column name. 
  cem.data <- as.data.frame(cbind(treat,X))
  
  mat <-
    cem(treatment="treat",data=cem.data,verbose=as.integer(verbose)+1,
        method=k2k.method,...)

  # here we create a column vector where the matched entry get its stratum
  # and the unmatched entry gets an NA.
  strat <- rep(NA,n)
  names(strat) <- names(treat)
  strat[mat$matched] <- mat$strata[mat$matched]

  # here we just add the names onto the wieght from the cem output
  wh <- mat$w
  names(wh) <- names(treat)

  # weighting functions in matchit error-out on these conditions,
  # so we should too.
 
  if (sum(wh)==0) 
    stop("No units were matched")
  else if (sum(wh[treat==1])==0)
    stop("No treated units were matched")
  else if (sum(wh[treat==0])==0)
    stop("No control units were matched")
  
  res <- list(subclass = strat, weights = mat$w)
  class(res) <- "matchit"
  return(res)
}
