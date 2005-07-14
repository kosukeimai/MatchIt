# Need to account for weights -- how do we do qq plots with weights
plot.matchit <- function(x,discrete.cutoff=5,type="QQ"){
  if (class(x)[1]=="matchit.exact"){
    return(warning("Plot() not appropriate for exact matching.  No plots generated"))
  }
  if(type=="QQ"){
    matchit.qqplot(x,discrete.cutoff)
  } else if(type=="jitter"){
    jitter.pscore(x)
  } else {
    stop("Invalid type")
  }
}
