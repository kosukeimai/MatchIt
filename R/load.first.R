.onAttach <- function(...) {
  mylib <- dirname(system.file(package = "MatchIt"))
  ver <- packageDescription("MatchIt", lib = mylib)$Version
  builddate <- packageDescription("MatchIt", lib = mylib)$Date
  cat(paste("## \n##  MatchIt (Version ", ver, ", built: ", builddate, ")\n", sep = "")) 
  cat("##  Please refer to http://gking.harvard.edu/matchit for full documentation \n",
      "##  or help.matchit() for help with commands supported by MatchIt.\n##\n",
      sep="")
}
