match.data <- function(object, group = "all") {
  data <- object$data
  treat <- eval(object$treat, data)
  if (group == "all")
    return(data[object$matched,])
  else if (group == "treat")
    return(data[object$matched & treat == 1,])
  else if (group == "control")
    return(data[object$matched & treat == 0,])
  else
    stop("error: invalid input for group.")
}
