matchit <- function(formula, data, method = "nearest", distance = "logit",
                    distance.options=list(), discard = "none",
                    reestimate = FALSE, ...) { 

  #Checking input format
  #data input
  if(is.null(data)) stop("Dataframe must be specified",call.=FALSE)
  if(!is.data.frame(data)){
    stop("Data must be a dataframe",call.=FALSE)}
  if(sum(is.na(data))>0)
    stop("Missing values exist in the data")

  ## 7/13/06: Convert character variables to factors as necessary
  ischar <- rep(0, dim(data)[2])
  for (i in 1:dim(data)[2]) 
	if(is.character(data[,i])) data[,i] <- as.factor(data[,i])
  
  ## check inputs
  if (!is.numeric(distance)) {
    fn1 <- paste("distance2", distance, sep = "")
    if (!exists(fn1))
      stop(distance, "not supported.")
  }
  fn2 <- paste("matchit2", method, sep = "")
  if (!exists(fn2))
    stop(method, "not supported.")

  ## obtain T and X
  tt <- terms(formula)
  attr(tt, "intercept") <- 0
  mf <- model.frame(tt, data)
  treat <- model.response(mf)
  X <- model.matrix(tt, data=mf)

  ## estimate the distance measure
  if (method == "exact") {
    distance <- out1 <- discarded <- NULL
    if (!is.null(distance))
      warning("distance is set to `NULL' when exact matching is used.")
  } 
  else if (is.numeric(distance)){
    out1 <- NULL
    discarded <- discard(treat, distance, discard, X)
  } else {
    if (is.null(distance.options$formula))
      distance.options$formula <- formula
    if (is.null(distance.options$data))
      distance.options$data <- data
    out1 <- do.call(fn1, distance.options)
    discarded <- discard(treat, out1$distance, discard, X)
    if (reestimate) {
      distance.options$data <- data[!discarded,]
      distance.options$weights <- distance.options$weights[!discard]
      tmp <- out1
      out1 <- do.call(fn1, distance.options)
      tmp$distance[!discarded] <- out1$distance
      out1$distance <- tmp$distance
    }
    distance <- out1$distance
  }

  ## matching!
  out2 <- do.call(fn2, list(treat, X, data, distance=distance, discarded, ...)) 
  
  ## putting all the results together
  out2$call <- match.call()
  out2$model <- out1$model
  out2$formula <- formula
  out2$treat <- treat
  if (is.null(out2$X)){
    out2$X <- X
  }
  out2$distance <- distance
  out2$discarded <- discarded
  
  return(out2)
}
