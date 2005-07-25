# Need to account for weights -- how do we do qq plots with weights
plot.matchit <- function(x, y = NULL, discrete.cutoff=5,type="QQ",numdraws=5000){
  if (class(x)[1]=="matchit.exact"){
    return(warning("Plot() not appropriate for exact matching.  No plots generated."))
  }
  if(type=="QQ"){
    matchit.qqplot(x=x,discrete.cutoff=discrete.cutoff,numdraws=numdraws)
  } else if(type=="jitter"){
    jitter.pscore(x)
  } else {
    stop("Invalid type")
  }
}
