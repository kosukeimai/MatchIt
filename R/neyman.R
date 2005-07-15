neyman<-function(Y, object, bootstrap=NULL, counter=TRUE){

  mf <- match.call()
  data<-eval(object$call$data, sys.parent())
  Y <- eval(mf$Y, data)
  
  objclass <- class(object)
  if(objclass[1]=="matchit.exact"){
    if(!is.null(bootstrap)){
      warning("Bootstrap not used for exact matching")
    }
    res <- neyman.matchit.exact(Y, object)
  } else if("matchit"%in%objclass){
    res <- neyman.matchit(Y, object, bootstrap, counter)
  } else {
    stop("Input object is not a matchit object")
  }
  res
}
