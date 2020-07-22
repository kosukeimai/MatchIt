#' @export
match.data <- function(object, group = "all", distance = "distance",
                       weights = "weights", subclass = "subclass",
                       data = NULL, drop.unmatched = TRUE) {

  if (is.null(data)) {
    if (!is.null(object$model)) {
      env <- attributes(terms(object$model))$.Environment
    } else {
      env <- parent.frame()
    }
    data <- eval(object$call$data, envir = env)
  }
  else {
    if (!is.data.frame(data)) {
      if (is.matrix(data)) data <- as.data.frame.matrix(data)
      else stop("data must be a data frame.", call. = FALSE)
    }
    if (nrow(data) != length(object$treat)) {
      stop("data must have as many rows as there were units in the original call to matchit().", call. = FALSE)
    }
  }

  vars <- names(data)

  if (!is.null(object$distance)) {
    if (distance %in% vars) {
      if (distance == "distance") {
        stop(paste0("\"", distance, "\" is already the name of a variable in the data. Please choose another name for distance using the 'distance' argument."), call. = FALSE)
      }
      else {
        stop(paste0("\"", distance, "\" is already the name of a variable in the data. Please choose another name for distance."), call. = FALSE)
      }
    }

    data[[distance]] <- object$distance
  }

  if (!is.null(object$weights)){
    if (weights %in% vars) {
      if (weights == "weights") {
        stop(paste0("\"", weights, "\" is already the name of a variable in the data. Please choose another name for weights using the 'weights' argument."), call. = FALSE)
      }
      else {
        stop(paste0("\"", weights, "\" is already the name of a variable in the data. Please choose another name for weights."), call. = FALSE)
      }
    }
    data[[weights]] <- object$weights
  }

  if (!is.null(object$subclass)){
    if (subclass %in% vars) {
      if (subclass == "subclass") {
        stop(paste0("\"", subclass, "\" is already the name of a variable in the data. Please choose another name for subclass using the 'subclass' argument."), call. = FALSE)
      }
      else {
        stop(paste0("\"", subclass, "\" is already the name of a variable in the data. Please choose another name for subclass."), call. = FALSE)
      }
    }
    data[[subclass]] <- object$subclass
  }

  treat <- object$treat
  if (drop.unmatched && !is.null(object$weights)) {
    data <- data[object$weights > 0,]
    treat <- treat[object$weights > 0]
  }

  group <- match_arg(group, c("all", "treated", "control"))
  if (group == "all") return(data)
  else if (group == "treated") return(data[treat == 1,])
  else if (group == "control") return(data[treat == 0,])
}
