neyman<-function(Y, object, bootstrap = NULL, verbose = TRUE){

  mf <- match.call()
  data<-eval(object$call$data, sys.parent())
  Y <- eval(mf$Y, data)
  
  objclass <- class(object)
  if(objclass[1]%in%c("matchit.exact",
                      "matchit.full")){
    if(!is.null(bootstrap)){
      warning("Bootstrap not used for exact or full matching")
    }
    res <- neyman.matchit.exactfull(Y, object)
  } else if("matchit"%in%objclass){
    res <- neyman.matchit(Y, object, bootstrap, verbose)
  } else {
    stop("Input object is not a matchit object")
  }
  res
}
