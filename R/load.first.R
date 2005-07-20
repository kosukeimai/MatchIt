.onAttach <- function(...) {
  cat("\nPlease refer to http://gking.harvard.edu/matchit for full documentation \n",
      "or help.matchit() for help with commands and matching methods supported by MatchIt.\n\n",
      sep='')
  if(!any(search()=="package:MASS"))
    require(MASS) 
}
