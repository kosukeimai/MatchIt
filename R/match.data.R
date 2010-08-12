match.data <- function(object, group = "all", distance = "distance",
                       weights = "weights", subclass = "subclass") {

  if (!is.null(object$model)) {
    env <- attributes(terms(object$model))$.Environment
  } else {
    env <- parent.frame()
  }
  data <- eval(object$call$data, envir = env)
  treat <- object$treat
  wt <- object$weights
  vars <- names(data)
  if (distance %in% vars)
    stop("invalid input for distance. choose a different name.")
  else if (!is.null(object$distance)) {
    dta <- data.frame(cbind(data, object$distance))
    names(dta) <- c(names(data), distance)
    data <- dta
  }
  if (weights %in% vars)
    stop("invalid input for weights. choose a different name.")
  else if (!is.null(object$weights)){
    dta <- data.frame(cbind(data, object$weights))
    names(dta) <- c(names(data), weights)
    data <- dta
  }
  if (subclass %in% vars)
    stop("invalid input for subclass. choose a different name.")
  else if (!is.null(object$subclass)){
    dta <- data.frame(cbind(data, object$subclass))
    names(dta) <- c(names(data), subclass)
    data <- dta
  }
  if (group == "all")
    return(data[wt > 0,])
  else if (group == "treat")
    return(data[wt > 0 & treat == 1,])
  else if (group == "control")
    return(data[wt > 0 & treat == 0,])
  else
    stop("error: invalid input for group.")
}
