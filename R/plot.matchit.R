# Need to account for weights -- how do we do qq plots with weights
plot.matchit <- function(x, discrete.cutoff=5, type="QQ",
                         numdraws=5000, interactive = T, ...){
  if ("matchit.exact" %in% class(x)){
    stop("Not appropriate for exact matching.  No plots generated.")
  }
  if(type=="QQ"){
    matchit.qqplot(x=x,discrete.cutoff=discrete.cutoff,
                   numdraws=numdraws,interactive=interactive)
  } else if(type=="jitter"){
    jitter.pscore(x)
  } else {
    stop("Invalid type")
  }
}
