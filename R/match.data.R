match.data <- function(object, group = "all", distance = "distance", weights = "weights", subclass = "subclass",
                               data = NULL, include.s.weights = TRUE, drop.unmatched = TRUE) {

  if (!inherits(object, "matchit")) {
    stop("'object' must be a matchit object, the output of a call to matchit().", call. = FALSE)
  }

  if (is.null(data)) {
    env <- environment(object$formula)
    data <- try(eval(object$call$data, envir = env), silent = TRUE)
    if (length(data) == 0 || inherits(data, "try-error") || length(dim(data)) != 2 || nrow(data) != length(object[["treat"]])) {
      env <- parent.frame()
      data <- try(eval(object$call$data, envir = env), silent = TRUE)
      if (length(data) == 0 || inherits(data, "try-error") || length(dim(data)) != 2 || nrow(data) != length(object[["treat"]])) {
        data <- object[["model"]][["data"]]
        if (length(data) == 0 || nrow(data) != length(object[["treat"]])) {
          stop("A valid dataset could not be found. Please supply an argument to 'data' containing the original dataset used in the matching.", call. = FALSE)
        }
      }
    }
  }

  if (!is.data.frame(data)) {
    if (is.matrix(data)) data <- as.data.frame.matrix(data)
    else stop("'data' must be a data frame.", call. = FALSE)
  }
  if (nrow(data) != length(object$treat)) {
    stop("'data' must have as many rows as there were units in the original call to matchit().", call. = FALSE)
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

  if (!is.null(object$weights)) {
    if (is.null(weights)) stop("The argument to 'weights' cannot be NULL.", call. = FALSE)
    if (!is.atomic(weights) || !is.character(weights) || length(weights) != 1 || is.na(weights)) {
      stop("The argument to 'weights' must be a string of length 1.", call. = FALSE)
    }
    if (weights %in% names(data)) {
      stop(paste0("\"", weights, "\" is already the name of a variable in the data. Please choose another name for weights using the 'weights' argument."), call. = FALSE)
    }
    data[[weights]] <- object$weights

    if (!is.null(object$s.weights) && include.s.weights) {
      data[[weights]] <- data[[weights]] * object$s.weights
    }
  }

  if (!is.null(object$subclass)) {
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
    data <- data[object$weights > 0,,drop = FALSE]
    treat <- treat[object$weights > 0]
  }

  group <- match_arg(group, c("all", "treated", "control"))
  if (group == "treated") data <- data[treat == 1,,drop = FALSE]
  else if (group == "control") data <- data[treat == 0,]

  if (!is.null(object$distance)) attr(data, "distance") <- distance
  if (!is.null(object$weights)) attr(data, "weights") <- weights
  if (!is.null(object$subclass)) attr(data, "subclass") <- subclass

  class(data) <- c("matchdata", class(data))

  return(data)
}

get_matches <- function(object, distance = "distance", weights = "weights", subclass = "subclass",
                                id = "id", data = NULL, include.s.weights = TRUE) {

  if (!inherits(object, "matchit")) {
    stop("'object' must be a matchit object, the output of a call to matchit().", call. = FALSE)
  }

  if (is.null(object$match.matrix)) {
    stop("A match.matrix component must be present in the matchit object, which does not occur with all types of matching. Use match.data() instead.", call. = FALSE)
  }

  #Get initial data using match.data; note weights and subclass will be removed,
  #including them here just checks their names don't clash
  m.data <- match.data(object, group = "all", distance = distance,
                       weights = weights, subclass = subclass, data = data,
                       include.s.weights = FALSE, drop.unmatched = TRUE)

  if (is.null(id)) stop("The argument to 'id' cannot be NULL.", call. = FALSE)
  if (!is.atomic(id) || !is.character(id) || length(id) != 1 || is.na(id)) {
    stop("The argument to 'id' must be a string of length 1.", call. = FALSE)
  }

  if (id %in% names(m.data)) {
    stop(paste0("\"", id, "\" is already the name of a variable in the data. Please choose another name for id using the 'id' argument."), call. = FALSE)
  }

  m.data[[id]] <- names(object$treat)[object$weights > 0]

  for (i in c(weights, subclass)) {
    if (i %in% names(m.data)) m.data[[i]] <- NULL
  }

  mm <- object$match.matrix
  mm <- mm[!is.na(mm[,1]),,drop = FALSE]
  tmm <- t(mm)

  num.matches <- rowSums(!is.na(mm))

  matched <- as.data.frame(matrix(NA_character_, nrow = nrow(mm) + sum(!is.na(mm)), ncol = 3))
  names(matched) <- c(id, subclass, weights)

  matched[[id]] <- c(as.vector(tmm[!is.na(tmm)]), rownames(mm))
  matched[[subclass]] <- c(as.vector(col(tmm)[!is.na(tmm)]), seq_len(nrow(mm)))
  matched[[weights]] <- c(1/num.matches[matched[[subclass]][seq_len(sum(!is.na(mm)))]], rep(1, nrow(mm)))

  if (!is.null(object$s.weights) && include.s.weights) {
    matched[[weights]] <- matched[[weights]] * object$s.weights[matched[[id]]]
  }

  out <- merge(matched, m.data, by = id, all.x = TRUE, sort = FALSE)

  out <- out[order(out[[subclass]], object$treat[out[[id]]], method = "radix", decreasing = c(FALSE, TRUE)),]
  rownames(out) <- NULL

  out[[subclass]] <- factor(out[[subclass]], labels = seq_len(nrow(mm)))

  if (!is.null(object$distance)) attr(out, "distance") <- distance
  attr(out, "weights") <- weights
  attr(out, "subclass") <- subclass
  attr(out, "id") <- id

  class(out) <- c("getmatches", class(out))

  return(out)
}

rbind.matchdata <- function(..., deparse.level = 1) {

  allargs <- list(...)
  allargs <- allargs[lengths(allargs) > 0L]
  if (is.null(names(allargs))) {
    md_list <- allargs
    allargs <- list()
  }
  else {
    md_list <- allargs[names(allargs) == ""]
    allargs[names(allargs) == ""] <- NULL
  }
  allargs$deparse.level <- deparse.level

  type <- intersect(c("matchdata", "getmatches"), unlist(lapply(md_list, class)))
  if (length(type) == 0) stop("A matchdata or getmatches object must be supplied.", call. = FALSE)
  else if (length(type) == 2) stop("Supplied objects must be all matchdata objects or all getmatches objects.", call. = FALSE)

  attrs <- c("distance", "weights", "subclass", "id")
  attr_list <- setNames(vector("list", length(attrs)), attrs)
  key_attrs <- setNames(rep(NA_character_, length(attrs)), attrs)

  for (i in attrs) {
    attr_list[[i]] <- unlist(lapply(md_list, function(m) {
      a <- attr(m, i)
      if (length(a) == 0) NA_character_ else a
    }))
    if (all(is.na(attr_list[[i]]))) attr_list[[i]] <- NULL
    else key_attrs[i] <- attr_list[[i]][which(!is.na(attr_list[[i]]))[1]]
  }
  attrs <- names(attr_list)
  key_attrs <- key_attrs[attrs]

  #Check if all non-attr columns are the same across datasets
  other_col_list <- lapply(seq_along(md_list), function(d) {
    setdiff(names(md_list[[d]]), unlist(lapply(attr_list, `[`, d)))
  })
  for (d in seq_along(md_list)[-1]) {
    if (length(other_col_list[[d]]) != length(other_col_list[[1]]) || !all(other_col_list[[d]] %in% other_col_list[[1]])) {
      stop(paste("The", switch(type, "matchdata" = "match.data()", "get_matches()"), "inputs must come from the same dataset."), call. = FALSE)
    }
  }

  for (d in seq_along(md_list)) {
    for (i in attrs) {
      #Rename columns of each attribute the same across datasets
      if (is.null(attr(md_list[[d]], i))) md_list[[d]] <- setNames(cbind(md_list[[d]], NA), c(names(md_list[[d]]), key_attrs[i]))
      else names(md_list[[d]])[names(md_list[[d]]) == attr_list[[i]][d]] <- key_attrs[i]

      #Give subclasses unique values across datasets
      if (i == "subclass") {
        if (all(is.na(md_list[[d]][[key_attrs[i]]]))) md_list[[d]][[key_attrs[i]]] <- factor(md_list[[d]][[key_attrs[i]]], levels = NA)
        else levels(md_list[[d]][[key_attrs[i]]]) <- paste(d, levels(md_list[[d]][[key_attrs[i]]]), sep = "_")
      }
    }

    #Put all columns in the same order
    if (d > 1) {
      md_list[[d]] <- md_list[[d]][names(md_list[[1]])]
    }
    class(md_list[[d]]) <- class(md_list[[d]])[class(md_list[[d]]) != type]
  }



  out <- do.call("rbind", c(md_list, allargs))

  for (i in attrs) {
    attr(out, i) <- unname(key_attrs[i])
  }

  class(out) <- c(type, class(out))

  out
}
rbind.getmatches <- rbind.matchdata