# Need to account for weights -- how do we do qq plots with weights
plot.matchit <- function(obj,discrete.cutoff=5,type="QQ"){
  if (class(obj)[1]=="matchit.exact"){
    return(warning("Plot() not appropriate for exact matching.  No plots generated"))
  }
  if(type=="QQ"){
    matchit.qqplot(obj,discrete.cutoff)
  } else if(type=="jitter"){
    jitter.pscore(obj)
  }
}
