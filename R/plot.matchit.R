# Need to account for weights -- how do we do qq plots with weights
plot.matchit <- function(x, discrete.cutoff=5, type="QQ",
                         numdraws=5000, interactive = T, which.xs =
                         NULL, ...){
  if ("matchit.exact" %in% class(x)){
    stop("Not appropriate for exact matching.  No plots generated.")
  }
  if(type=="QQ"){
    matchit.qqplot(x=x,discrete.cutoff=discrete.cutoff,
                   numdraws=numdraws, interactive=interactive,
                   which.xs = which.xs, ...)
  } else if(type=="jitter"){
    if("matchit.mahalanobis" %in% class(x)){
      stop("Not appropriate for pure Mahalanobis matching.  No plots generated.")
    }
    jitter.pscore(x, interactive=interactive,...)
  } else if(type=="hist"){
    if("matchit.mahalanobis" %in% class(x)){
      stop("Not appropriate for pure Mahalanobis matching.  No plots generated.")
    }
    hist.pscore(x,...)
  } else {
    stop("Invalid type")
  }
}
