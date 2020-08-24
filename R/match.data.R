#' @export
match.data <- function(object, group = "all", distance = "distance",
                       weights = "weights", subclass = "subclass",
                       data = NULL, drop.unmatched = TRUE) {

  if (!inherits(object, "matchit")) {
    stop("'object' must be a matchit object, the output of a call to matchit().", call. = FALSE)
  }
  if (is.null(data)) {
    if (!is.null(object$model)) {
      env <- attributes(terms(object$model))$.Environment
    } else {
      env <- parent.frame()
    }
    data <- eval(object$call$data, envir = env)
    if (length(data) == 0) stop("A dataset could not be found. Please supply an argument to 'data' containing the original dataset used in the matching.", call. = FALSE)
  }
  else {
    if (!is.data.frame(data)) {
      if (is.matrix(data)) data <- as.data.frame.matrix(data)
      else stop("'data' must be a data frame.", call. = FALSE)
    }
    if (nrow(data) != length(object$treat)) {
      stop("'data' must have as many rows as there were units in the original call to matchit().", call. = FALSE)
    }
  }

  if (!is.null(object$distance)) {
    if (is.null(distance)) stop("The argument to 'distance' cannot be NULL.", call. = FALSE)
    if (!is.atomic(distance) || !is.character(distance) || length(distance) != 1 || is.na(distance)) {
      stop("The argument to 'distance' must be a string of length 1.", call. = FALSE)
    }
    if (distance %in% names(data)) {
      stop(paste0("\"", distance, "\" is already the name of a variable in the data. Please choose another name for distance using the 'distance' argument."), call. = FALSE)
    }
    data[[distance]] <- object$distance
  }

  if (!is.null(object$weights)){
    if (is.null(weights)) stop("The argument to 'weights' cannot be NULL.", call. = FALSE)
    if (!is.atomic(weights) || !is.character(weights) || length(weights) != 1 || is.na(weights)) {
      stop("The argument to 'weights' must be a string of length 1.", call. = FALSE)
    }
    if (weights %in% names(data)) {
      stop(paste0("\"", weights, "\" is already the name of a variable in the data. Please choose another name for weights using the 'weights' argument."), call. = FALSE)
    }
    data[[weights]] <- object$weights
  }

  if (!is.null(object$subclass)){
    if (is.null(subclass)) stop("The argument to 'subclass' cannot be NULL.", call. = FALSE)
    if (!is.atomic(subclass) || !is.character(subclass) || length(subclass) != 1 || is.na(subclass)) {
      stop("The argument to 'subclass' must be a string of length 1.", call. = FALSE)
    }
    if (subclass %in% names(data)) {
      stop(paste0("\"", subclass, "\" is already the name of a variable in the data. Please choose another name for subclass using the 'subclass' argument."), call. = FALSE)
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

#' @export
get_matches <- function(object, distance = "distance", weights = "weights",
                        subclass = "subclass", id = "id", data = NULL) {

  if (!inherits(object, "matchit")) {
    stop("'object' must be a matchit object, the output of a call to matchit().", call. = FALSE)
  }
  if (is.null(object$match.matrix)) {
    stop("A match.matrix component must be present in the matchit object, which does not occur with all types of matching. Use match.data() instead.", call. = FALSE)
  }

  m.data <- match.data(object, group = "all", distance = distance,
                       weights = weights, subclass = subclass, data = data,
                       drop.unmatched = TRUE)

  if (is.null(id)) stop("The argument to 'id' cannot be NULL.", call. = FALSE)
  if (!is.atomic(id) || !is.character(id) || length(id) != 1 || is.na(id)) {
    stop("The argument to 'id' must be a string of length 1.", call. = FALSE)
  }
  if (id %in% names(data)) {
    stop(paste0("\"", id, "\" is already the name of a variable in the data. Please choose another name for id using the 'id' argument."), call. = FALSE)
  }

  m.data[[id]] <- rownames(m.data)

  m.data[c(weights, subclass)] <- NULL

  mm <- object$match.matrix

  mm <- mm[!is.na(mm[,1]),,drop = FALSE]

  tmm <- t(mm)

  num.matches <- rowSums(!is.na(mm))

  # #Create match matrix containing just the treated units (in the rownames)
  # mm1 <- matrix(rownames(mm), nrow = nrow(mm), ncol = ncol(mm))
  # mm1[is.na(mm)] <- NA_character_
  # tmm1 <- t(mm1)
  #
  #
  # matched <- as.data.frame(matrix(NA_character_, nrow = 2*sum(!is.na(mm)), ncol = 3))
  # names(matched) <- c(id, subclass, weights)
  #
  # matched[[id]] <- c(as.vector(tmm[!is.na(tmm)]), as.vector(tmm1[!is.na(tmm1)]))
  # matched[[subclass]] <- c(as.vector(col(tmm)[!is.na(tmm)]), as.vector(col(tmm1)[!is.na(tmm1)]))
  # matched[[weights]] <- c(1/num.matches[matched$subclass])

  matched <- as.data.frame(matrix(NA_character_, nrow = nrow(mm) + sum(!is.na(mm)), ncol = 3))
  names(matched) <- c(id, subclass, weights)

  matched[[id]] <- c(as.vector(tmm[!is.na(tmm)]), rownames(mm))
  matched[[subclass]] <- c(as.vector(col(tmm)[!is.na(tmm)]), seq_len(nrow(mm)))
  matched[[weights]] <- c(1/num.matches[matched$subclass[seq_len(sum(!is.na(mm)))]], rep(1, nrow(mm)))

  out <- merge(matched, m.data, by = id, all.x = TRUE, sort = FALSE)

  out <- out[order(out[[subclass]], object$treat[out[[id]]], method = "radix", decreasing = c(FALSE, TRUE)),]
  rownames(out) <- NULL

  return(out)
}