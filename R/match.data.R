match.data <- function(object, group = "all", distance = "distance",
                       weights = "weights") {
  data <- eval(object$call$data)
  treat <- object$treat
  vars <- names(data)
  if (distance %in% vars)
    stop("invalid input for distance. choose a different name.")
  if (weights %in% vars)
    stop("invalid input for weights. choose a different name.")
  if (group == "all")
    return(data[object$matched,])
  else if (group == "treat")
    return(data[object$matched & treat == 1,])
  else if (group == "control")
    return(data[object$matched & treat == 0,])
  else
    stop("error: invalid input for group.")
}
