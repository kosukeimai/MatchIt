
#' @import MASS
#' @import stats
#' @import utils
#' @import graphics

in_frame <- function(df, col) {
  exists(col, envir= as.environment(df))
}

int_rownames <- function(mat) {
  options(warn= -1)
  rn <- as.integer(rownames(mat))
  options(warn= 0)
  if (any(is.na(rn))) return(FALSE)
  else return(TRUE)
}

#----------------------------------------------------------
### Inheritance
#----------------------------------------------------------

#' @title Checks matchit Class
#' @description Function that checks if the target object is a \code{matchit} object.
#' @param object any R object
#' @return Returns \code{TRUE} if its argument has class "matchit" among its classes and
#' \code{FALSE} otherwise.
#' @export
is.matchit <- function(object) {
  inherits(object, "matchit")
}


#----------------------------------------------------------
### GET MATCHES
#----------------------------------------------------------

#' @title Get matches from matchit object
#' @description Get the resulting matches from a \code{matchit} model object. This function allows the
#' user to extract the matches from the original dataset used in model building or from a new dataset
#' that has a matching set of key column(s) (\code{id_cols}).
#' @param object The \code{'matchit'} class model object
#' @param model_frame The \code{'data.frame'} class object used in creation of \code{object}.
#' @param id_cols A string indicating the ID for the datset used in the call to \code{\link{matchit}}.
#' This can be used in combination with \code{newdata} to return the base dataset. Defaults to
#' \code{NULL}.
#' @param newdata A new \code{data.frame} object to extract matched observations from. Used in
#' conjunction with \code{id_cols}. Defaults to \code{NULL}.
#' @return If \code{newdata} is \code{NULL}, a subset of \code{model_frame} containing the rows
#' corresponding to the matched treatment and control observations with weights appended. If
#' \code{newdata} is not \code{NULL}, an equivalent subset of \code{newdata} is returned.
#' @export
get_matches <- function(object, model_frame, id_cols= NULL, newdata= NULL) {
  UseMethod("get_matches", object)
}


#' @export
get_matches.matchit <- function(object, model_frame, id_cols= NULL, newdata= NULL) {
  # 00. error checking
  if (!is.matchit(object)) stop("object must be a class 'matchit' object.")
  if (!is.data.frame(model_frame)) stop("model_frame must be a data.frame")
  if (nrow(model_frame) != nrow(object$X))
    stop("model_frame must have the same number of rows as the data used in the call to 'matchit'.")
  if ( (!is.null(id_cols) & is.null(newdata)) | (is.null(id_cols) & !is.null(newdata)) ) {
    stop("For identity returns, both id_cols and newdata must be supplied.")
  }
  if (!is.null(id_cols) & !is.null(newdata)) {
    if (!is.character(id_cols)) stop("id_cols must be a character vector.")
    if (!all(sapply(id_cols, in_frame, df= model_frame)))
      stop("all values in id_cols must exist in model_frame")
    if (!all(sapply(id_cols, in_frame, df= newdata)))
      stop("all values in id_cols must exist in newdata")
  }

  # 01. preliminaries
  if ( any(grepl(pattern= "replace", x= names(object$call))) ) {
    idx_call_replace <- which(names(object$call) == "replace")
  } else {
    idx_call_replace <- -1L
  }
  use_subclass_matching <- base::grepl(x = object$call[4],
                             pattern= paste(c("exact", "full", "subclass", "cem"), collapse= "|"),
                             ignore.case= TRUE)
  # [deprecated] use_genetic_matching <- base::grepl(x = object$call[4], pattern= "genetic", ignore.case= TRUE)
  use_newdata <- ifelse(is.null(newdata), FALSE, TRUE)
  has_int_rownames <- int_rownames(model_frame)

  # 02. extract pairs, either exact matching or otherwise
  if (use_subclass_matching) {
    if (has_int_rownames) {
      row_idx <- as.integer(names(object$subclass)[which(!is.na(object$subclass))])
      match_wts <- object$weights[which(!is.na(object$subclass))]
      model_subset <- data.frame(model_frame[row_idx, ], weight= match_wts)
    } else {
      row_idx <- names(object$subclass)[which(!is.na(object$subclass))]
      match_wts <- object$weights[which(!is.na(object$subclass))]
      model_subset <- data.frame(model_frame[which(rownames(model_frame) %in% row_idx), ],
                                 weight= match_wts)
    }
  } else {
   model_subset <- get_matches_non_subclass(object= object, model_frame= model_frame,
                     has_int_rownames= has_int_rownames,
                     idx_call_replace= idx_call_replace)
  }

  # 03. return
  if (!use_newdata) {
    return(model_subset)
  } else { # using newdata via id_cols
    newdata <- base::as.data.frame(newdata) # in case of data.table
    unique_ids <- base::unique(x= model_subset[, c(id_cols, "weight")])
    # use of all == FALSE for two cases:
    # a) newdata contains additional IDs not in the matching data -- most common
    # b) if model_subset contains IDs that newdata does not have -- less common
    return(merge(newdata, unique_ids, by= id_cols, all.x=FALSE, all.y=FALSE))
  }
}

get_matches_non_subclass <- function(object, model_frame, has_int_rownames, use_genetic_matching,
                                     idx_call_replace) {
  ## get match weights
  n_matches <- ncol(object$match.matrix)
  # get control observations that have matches and their frequency (ie weight)
  control_units <- object$match.matrix
  attr(control_units, "dim") <- NULL
  control_units <- control_units[!is.na(control_units)]

  control_wts <- table(control_units); class(control_wts) <- "vector"
  control_wts <- data.frame(ob= names(control_wts), weight= control_wts / n_matches,
                            stringsAsFactors= FALSE)

  if (has_int_rownames) {
    treated_obs <- data.frame(model_frame[as.integer(rownames(object$match.matrix)), ], weight= 1)
    control_obs <- data.frame(model_frame[as.integer(control_wts$ob), ],
                              weight= control_wts$weight)
  } else {
    treated_obs <- data.frame(
      model_frame[which(rownames(model_frame) %in% rownames(object$match.matrix)), ],
      weight= 1)

    model_frame$rownames <- rownames(model_frame)
    control_obs <- merge(model_frame, control_wts, by.x= "rownames", by.y= "ob", all= FALSE)
    rownames(control_obs) <- control_obs$rownames
    control_obs$rownames <- NULL
    control_obs$ob <- NULL
  }

  model_subset <- do.call("rbind", list(control_obs, treated_obs))
  return(model_subset)
}

#----------------------------------------------------------
### PLOT METHODS
#----------------------------------------------------------

# Need to account for weights -- how do we do qq plots with weights
#' @export
plot.matchit <- function(x, discrete.cutoff = 5, type = "qq",
                         interactive = TRUE, which.xs = NULL, ...) {

  type <- tolower(type)
  type <- match_arg(type, c("qq", "jitter", "histogram"))

  if (type == "qq") {
    matchit.qqplot(x, discrete.cutoff=discrete.cutoff,
                   interactive=interactive,
                   which.xs = which.xs, ...)
  }
  else if (type == "jitter") {
    if (is.null(x$distance)) {
      stop("type = \"jitter\" cannot be used if a distance measure is not estimated or supplied. No plots generated.", call. = FALSE)
    }
    jitter.pscore(x, interactive = interactive,...)
  }
  else if (type=="histogram"){
    if (is.null(x$distance)){
      stop("type = \"hist\" cannot be used if a distance measure is not estimated or supplied. No plots generated.", call. = FALSE)
    }
    hist.pscore(x,...)
  }
}

#' @export
plot.matchit.subclass <- function(x, discrete.cutoff=5,
                                  type="qq", interactive = T,
                                  subclass = NULL, which.xs=NULL,...){
  choice.menu <- function(choices,question)
  {
    k <- length(choices)-1
    Choices <- data.frame(choices)
    row.names(Choices) <- 0:k
    names(Choices) <- "Choices"
    print.data.frame(Choices,right=FALSE)
    ans <- readline(question)
    while(!ans %in% 0:k) {
      print("Not valid -- please pick one of the choices")
      print.data.frame(Choices, right=FALSE)
      ans <- readline(question)
    }
    return(ans)
  }

  type <- tolower(type)
  type <- match_arg(type, c("qq", "jitter", "histogram"))

  if (type=="qq"){
    if (interactive) {
      subclasses <- sort(unique(x$subclass[!is.na(x$subclass)]))
      choices <- c("No", paste0("Yes : Subclass ", subclasses))
      question <- "Would you like to see quantile-quantile plots of any subclasses?"
      ans <- -1
      while(ans != 0) {
        ans <- as.numeric(choice.menu(choices, question))
        if (ans != 0) {
          matchit.qqplot(x,discrete.cutoff,which.subclass=subclasses[ans],
                         interactive = interactive, which.xs=which.xs,...)
        }
      }
    }
    else {
      matchit.qqplot(x, discrete.cutoff, which.subclass = subclass,
                     interactive = interactive, which.xs = which.xs,...)
    }
  }
  else if(type=="jitter"){
    jitter.pscore(x, interactive = interactive,...)
  }
  else if(type=="histogram"){
    hist.pscore(x,...)
  }
}


#' @export
plot.summary.matchit <- function(x, abs = TRUE, var.order = "data", threshold = c(.1, .05), position = "bottomright", ...) {
  if (!"Std. Mean Diff." %in% colnames(x$sum.all)) {
    stop("Not appropriate for unstandardized summary.  Run summary() with the standardize = TRUE option, and then plot.", call. = FALSE)
  }

  .pardefault <- par(no.readonly = TRUE)
  on.exit(par(.pardefault))

  sub <- inherits(x, "summary.matchit.subclass")

  var.order <- match_arg(var.order, c("data", "matched", "unmatched", "alphabetical"))

  sd.all <- x[["sum.all"]][,"Std. Mean Diff."]
  sd.matched <- x[[if (sub) "sum.across" else "sum.matched"]][,"Std. Mean Diff."]

  if (abs) {
    sd.all <- abs(sd.all)
    sd.matched <- abs(sd.matched)
    xlab <- "Absolute Standardized Mean Difference"
  }
  else {
    xlab <- "Standardized Mean Difference"
  }

  ord <- switch(var.order,
                "data" = rev(seq_along(sd.all)),
                "matched" = order(sd.matched),
                "unmatched" = order(sd.all),
                "alphabetical" = order(rownames(x$sum.all), decreasing = TRUE))

  dotchart(sd.all[ord], labels = rownames(x$sum.all)[ord], xlab = xlab,
           bg = NA, col = NA)
  abline(v = 0)

  if (sub && length(x$sum.subclass) > 0) {
    for (i in seq_along(x$sum.subclass)) {
      sd.sub <- x$sum.subclass[[i]][,"Std. Mean Diff."]
      if (abs) sd.sub <- abs(sd.sub)
      points(x = sd.sub[ord], y = seq_along(sd.sub),
             pch = as.character(i), col = "gray60", cex = .6)
    }
  }

  points(x = sd.all[ord], y = seq_along(sd.matched),
         pch = 21, bg = "white", col = "black")
  points(x = sd.matched[ord], y = seq_along(sd.matched),
         pch = 21, bg = "black", col = "black")

  if (!is.null(threshold)) {
    if (abs) {
        abline(v = threshold, lty = seq_along(threshold))
    }
    else {
      abline(v = threshold, lty = seq_along(threshold))
      abline(v = -threshold, lty = seq_along(threshold))
    }
  }

  if (!is.null(position)) {
    position <- match_arg(position, c("bottomright", "bottom", "bottomleft", "left",
                                      "topleft", "top", "topright", "right", "center"))
    legend(position, legend = c("All", "Matched"),
           pt.bg = c("white", "black"), pch = 21,
           inset = .015, xpd = T)
  }
  invisible(x)
}

#----------------------------------------------------------
### PRINT METHODS
#----------------------------------------------------------

#' @export
print.matchit <- function(x, digits = getOption("digits"), ...){
  cat("\nCall: ", deparse(x$call), sep="\n")
  cat("\nSample sizes:\n")

  #if(any(x$weights>0))
  #  nn <- rbind(table(x$treat),
  #              table(x$weights>0, x$treat),
  #              c(0,0))
  #else
  #  nn <- rbind(table(x$treat),
  #              table(x$weights>0,x$treat)[2:1,])

  print.table(x$nn, ...)
  invisible(x)
  cat("\n")
}

#' @export
print.matchit.exact <- function(x, digits = getOption("digits"), ...){
  cat("\nCall: ", deparse(x$call), sep = "\n")
  cat("\nExact Subclasses: ", length(unique(x$subclass[!is.na(x$subclass)])),"\n", sep="")
  cat("\nSample sizes:\n")
  # ntab <- table(factor(!is.na(x$subclass),
  #                      levels=c("TRUE","FALSE")),
  #               x$treat)
  # nn <- rbind(table(x$treat),
  #             ntab[c("TRUE","FALSE"),])
  # dimnames(nn) <- list(c("All","Matched","Unmatched"),
  #                      c("Control","Treated"))

  print.table(x$nn, ...)
  invisible(x)
  cat("\n")
}

#' @export
print.matchit.subclass <- function(x, digits = getOption("digits"), ...){
  cat("\nCall: ", paste(deparse(x$call), collapse = "\n"), sep = "\n")
  cat("\nSample sizes by subclasses:\n\n")

  qn <- table(x$treat[!x$discarded], x$subclass[!x$discarded])
  rownames(qn) <- c("Control", "Treated")

  if (any(x$discarded)) {
    qn <- cbind(qn, table(x$treat[x$discarded]))
    colnames(qn)[ncol(qn)] <- "Discarded"
  }
  qn <- rbind(qn, colSums(qn))
  rownames(qn)[nrow(qn)] <- "Total"

  qn <- cbind(qn, rowSums(qn))
  colnames(qn)[ncol(qn)] <- "All"

  print.table(qn, ...)
  invisible(x)
  cat("\n")
}

#' @export
print.summary.matchit <- function(x, digits = max(3, getOption("digits") - 3), ...){

  cat("\nCall:", deparse(x$call), sep = "\n")
  cat("\nSummary of balance for all data:\n")
  print.data.frame(as.data.frame(round(x$sum.all,digits)))
  cat("\n")

  if(!is.null(x$sum.matched)) {
    cat("\nSummary of balance for matched data:\n")
    print.data.frame(as.data.frame(round(x$sum.matched,digits)))
    cat("\nPercent Balance Improvement:\n")
    print.data.frame(as.data.frame(round(x$reduction,digits)))
    cat("\nSample sizes:\n")
    print.table(x$nn, digits=digits)
    cat("\n")
  }
  invisible(x)
}

#' @export
print.summary.matchit.subclass <- function(x, digits = max(3, getOption("digits") -  3), ...){
  cat("\nCall:", deparse(x$call), sep = "\n")
  cat("\nSummary of balance for all data:\n")
  print.data.frame(as.data.frame(round(x$sum.all,digits)))
  if (length(x$sum.subclass) > 0) {
    cat("\nSummary of balance by subclasses:\n")
    for (s in seq_along(x$sum.subclass)) {
      cat(paste0("\n- ", names(x$sum.subclass)[s], "\n"))
      print.data.frame(as.data.frame(round(x$sum.subclass[[s]], digits)))
    }
  }
  cat("\nSummary of balance across subclasses\n")
  print.data.frame(as.data.frame(round(x$sum.across, digits)))
  cat("\nPercent Balance Improvement:\n")
  print.data.frame(as.data.frame(round(x$reduction, digits)))
  cat("\nSample sizes by subclasses:\n")
  print.table(x$qn)
  cat("\n")
}


#----------------------------------------------------------
### SUMMARY METHODS
#----------------------------------------------------------

#' @export
summary.matchit <- function(object, interactions = FALSE,
                            addlvariables = NULL, standardize = FALSE,
                            data = NULL, ...) {

  #Create covariate matrix; include exact and mahvars
  X <- model.matrix(update(object$formula, NULL ~ . + 1), data = object$X,
                    contrasts.arg = lapply(Filter(is.factor, object$X),
                                           function(x) contrasts(x, contrasts = FALSE)))[,-1,drop = FALSE]
  if (!is.null(object$exact)) {
    Xexact <- model.matrix(update(object$exact, NULL ~ . + 1), data = object$X,
                           contrasts.arg = lapply(Filter(is.factor, object$X[, all.vars(object$exact), drop = FALSE]),
                                                  function(x) contrasts(x, contrasts = FALSE)))[,-1,drop = FALSE]
    X <- cbind(X, Xexact[,setdiff(colnames(Xexact), colnames(X)), drop = FALSE])
  }

  if (!is.null(object$mahvars)) {
    Xmahvars <- model.matrix(update(object$mahvars, NULL ~ . + 1), data = object$X,
                           contrasts.arg = lapply(Filter(is.factor, object$X[, all.vars(object$mahvars), drop = FALSE]),
                                                  function(x) contrasts(x, contrasts = FALSE)))[,-1,drop = FALSE]
    X <- cbind(X, Xmahvars[,setdiff(colnames(Xmahvars), colnames(X)), drop = FALSE])
  }

  if (!is.null(addlvariables)) {
    if (is.character(addlvariables)) {
      if (!is.null(data) && is.data.frame(data)) {
        if (all(addlvariables %in% names(data))) {
          addlvariables <- data[addlvariables]
        }
        else {
          stop("All variables in addlvariables must be in data.", call. = FALSE)
        }
      }
      else {
        stop("If addlvariables is specified as a string, a data frame argument must be supplied to data.", call. = FALSE)
      }
    }
    else if (inherits(addlvariables, "formula")) {
      vars.in.formula <- all.vars(addlvariables)
      if (!is.null(data) && is.data.frame(data)) data <- data.frame(data[names(data) %in% vars.in.formula], object$X[names(data) %in% setdiff(vars.in.formula, names(data))])
      else data <- object$X

      addlvariables <- model.matrix(update(addlvariables, NULL ~ . + 1), data = data,
                                    contrasts.arg = lapply(Filter(is.factor, data[, all.vars(addlvariables), drop = FALSE]),
                                                           function(x) contrasts(x, contrasts = FALSE)))[,-1,drop = FALSE]
    }
    else if (!is.matrix(addlvariables) && !is.data.frame(addlvariables)) {
      stop("The argument to addlvariables must be in one of the accepted forms. See ?summary.matchit for details.", call. = FALSE)
    }

    if (is.data.frame(addlvariables)) {
      if (!all(vapply(addlvariables, is.numeric, logical(1L)))) {
        addlvariables <- model.matrix(reformulate(names(addlvariables)), data = addlvariables,
                                      contrasts.arg = lapply(Filter(is.factor, addlvariables),
                                                             function(x) contrasts(x, contrasts = FALSE)))[,-1,drop = FALSE]

      }
      else {
        addlvariables <- as.matrix(addlvariables)
      }
    }

    X <- cbind(X, addlvariables[, setdiff(colnames(addlvariables), colnames(X)), drop = FALSE])
  }

  treat <- object$treat
  weights <- object$weights
  nam <- colnames(X)

  kk <- ncol(X)

  ## Summary Stats
  aa <- setNames(lapply(seq_len(kk), function(i) qoi(X[,i], tt = treat, ww = weights, standardize = standardize)),
                 colnames(X))

  sum.all <- sum.matched <- matrix(NA_real_, nrow = kk, ncol = ncol(aa[[1]]), dimnames = list(nam, colnames(aa[[1]])))

  sum.all.int <- sum.matched.int <- NULL
  for (i in 1:kk) {
    sum.all[i,] <- aa[[i]][1,]
    sum.matched[i,] <- aa[[i]][2,]
  }

  if (interactions) {
    sum.all.int <- sum.matched.int <- matrix(NA_real_, nrow = kk*(kk+1)/2, ncol = ncol(aa[[1]]), dimnames = list(NULL, colnames(aa[[1]])))
    to.remove <- rep(FALSE, nrow(sum.all.int))
    int.names <- character(nrow(sum.all.int))
    k <- 1
    for (i in 1:kk) {
      for (j in i:kk) {
        x2 <- X[,i] * X[,j]
        if (all(abs(x2) < sqrt(.Machine$double.eps)) ||
            all(abs(x2 - X[,i]) < sqrt(.Machine$double.eps))) { #prevent interactions within same factors
          to.remove[k] <- TRUE
        }
        else {
          jqoi <- qoi(x2, tt = treat, ww = weights, standardize = standardize)
          sum.all.int[k,] <- jqoi[1,]
          sum.matched.int[k,] <- jqoi[2,]
          if (i == j) {
            int.names[k] <- paste0(nam[i], "\u00B2")
          }
          else {
            int.names[k] <- paste(nam[i], nam[j], sep=" * ")
          }
        }
        k <- k + 1
      }
    }
    rownames(sum.all.int) <- rownames(sum.matched.int) <- int.names

    sum.all <- rbind(sum.all, sum.all.int[!to.remove,,drop = FALSE])
    sum.matched <- rbind(sum.matched, sum.matched.int[!to.remove,,drop = FALSE])
  }

  if (!is.null(object$distance)) {
    ad <- qoi(object$distance, tt = treat, ww = weights, standardize = standardize)
    sum.all <- rbind(ad[1,], sum.all)
    sum.matched <- rbind(ad[2,], sum.matched)
    rownames(sum.all)[1] <- rownames(sum.matched)[1] <- "distance"
  }

  ## Imbalance Reduction
  stat.all <- abs(sum.all[,-(1:2)])
  stat.matched <- abs(sum.matched[,-(1:2)])
  reduction <- 100*(stat.all-stat.matched)/stat.all

  reduction[stat.all == 0 & stat.matched == 0] <- 0
  reduction[stat.all == 0 & stat.matched > 0] <- -Inf

  ## output
  res <- list(call = object$call, nn = object$nn, sum.all = sum.all,
              sum.matched = sum.matched, reduction = reduction)
  class(res) <- "summary.matchit"
  return(res)
}

#' @export
summary.matchit.subclass <- function(object, interactions = FALSE,
                                     addlvariables = NULL, standardize = FALSE,
                                     data = NULL, subclass = FALSE, ...) {

  #Create covariate matrix; include exact and mahvars
  X <- model.matrix(update(object$formula, NULL ~ . + 1), data = object$X,
                    contrasts.arg = lapply(Filter(is.factor, object$X),
                                           function(x) contrasts(x, contrasts = FALSE)))[,-1,drop = FALSE]

  if (!is.null(addlvariables)) {
    if (is.character(addlvariables)) {
      if (!is.null(data) && is.data.frame(data)) {
        if (all(addlvariables %in% names(data))) {
          addlvariables <- data[addlvariables]
        }
        else {
          stop("All variables in addlvariables must be in data.", call. = FALSE)
        }
      }
      else {
        stop("If addlvariables is specified as a string, a data frame argument must be supplied to data.", call. = FALSE)
      }
    }
    else if (inherits(addlvariables, "formula")) {
      vars.in.formula <- all.vars(addlvariables)
      if (!is.null(data) && is.data.frame(data)) data <- data.frame(data[names(data) %in% vars.in.formula],
                                                                    object$X[names(data) %in% setdiff(vars.in.formula, names(data))])
      else data <- object$X

      addlvariables <- model.matrix(update(addlvariables, NULL ~ . + 1), data = data,
                                    contrasts.arg = lapply(Filter(is.factor, data[, all.vars(addlvariables), drop = FALSE]),
                                                           function(x) contrasts(x, contrasts = FALSE)))[,-1,drop = FALSE]
    }
    else if (!is.matrix(addlvariables) && !is.data.frame(addlvariables)) {
      stop("The argument to addlvariables must be in one of the accepted forms. See ?summary.matchit for details.", call. = FALSE)
    }

    if (is.data.frame(addlvariables)) {
      if (!all(vapply(addlvariables, is.numeric, logical(1L)))) {
        addlvariables <- model.matrix(reformulate(names(addlvariables)), data = addlvariables,
                                      contrasts.arg = lapply(Filter(is.factor, addlvariables),
                                                             function(x) contrasts(x, contrasts = FALSE)))[,-1,drop = FALSE]

      }
      else {
        addlvariables <- as.matrix(addlvariables)
      }
    }
  }

  X <- cbind(X, addlvariables[, setdiff(colnames(addlvariables), colnames(X)), drop = FALSE])

  which.subclass <- subclass
  treat <- object$treat
  weights <- object$weights
  subclass <- object$subclass
  nam <- colnames(X)

  kk <- ncol(X)
  subclasses <- sort(unique(subclass, nmax = length(object$q.cut) - 1))

  if (isTRUE(which.subclass)) which.subclass <- subclasses
  else if (isFALSE(which.subclass)) which.subclass <- NULL
  else if (!is.atomic(which.subclass) || !all(which.subclass %in% subclasses)) {
    stop("subclass should be TRUE, FALSE, or a vector of subclass indices for which subclass balance is to be displayed.")
  }
  else which.subclass <- subclasses[which.subclass]

  ## By Subclass
  sum.subclass <- lapply(which.subclass, function(s) {
    in.sub <- subclass == s & !object$discarded

    #qoi without weights only returns unmatched stats, which is all we need within
    #subclasses. Otherwise, identical to matched stats.
    aa <- setNames(lapply(seq_len(kk), function(i) qoi(X[in.sub,i], tt = treat[in.sub], standardize = standardize)),
             colnames(X))

    sum.sub <- matrix(NA_real_, nrow = kk, ncol = ncol(aa[[1]]), dimnames = list(nam, colnames(aa[[1]])))

    sum.sub.int <- NULL
    for (i in 1:kk) {
      sum.sub[i,] <- aa[[i]][1,]
    }
    if (interactions) {
      sum.sub.int <- matrix(NA_real_, nrow = kk*(kk+1)/2, ncol = ncol(aa[[1]]), dimnames = list(NULL, colnames(aa[[1]])))
      to.remove <- rep(FALSE, nrow(sum.sub.int))
      int.names <- character(nrow(sum.sub.int))
      k <- 1
      for (i in 1:kk) {
        for (j in i:kk) {
          x2 <- X[,i] * X[,j]
          if (all(abs(x2) < sqrt(.Machine$double.eps)) ||
              all(abs(x2 - X[,i]) < sqrt(.Machine$double.eps))) { #prevent interactions within same factors
            to.remove[k] <- TRUE
          }
          else {
            jqoi <- qoi(x2[in.sub], tt = treat[in.sub], standardize = standardize)
            sum.sub.int[k,] <- jqoi[1,]
            if (i == j) {
              int.names[k] <- paste0(nam[i], "\u00B2")
            }
            else {
              int.names[k] <- paste(nam[i], nam[j], sep=" * ")
            }
          }
          k <- k + 1
        }
      }
      rownames(sum.sub.int) <- int.names

      sum.sub <- rbind(sum.sub, sum.sub.int[!to.remove,,drop = FALSE])
    }

    if (!is.null(object$distance)) {
      ad <- qoi(object$distance[in.sub], tt = treat[in.sub], standardize = standardize)
      sum.sub <- rbind(ad[1,], sum.sub)
      rownames(sum.sub)[1] <- "distance"
    }

    return(sum.sub)
  })
  if (!is.null(which.subclass)) names(sum.subclass) <- paste("Subclass", which.subclass)

  ## Aggregate Subclass
  #Use the estimated weights to compute aggregate balance.
  ## Summary Stats
  aa <- setNames(lapply(seq_len(kk), function(i) qoi(X[,i], tt = treat, ww = weights, standardize = standardize)),
                 colnames(X))

  sum.all <- sum.matched <- matrix(NA_real_, nrow = kk, ncol = ncol(aa[[1]]), dimnames = list(nam, colnames(aa[[1]])))

  sum.all.int <- sum.matched.int <- NULL
  for (i in 1:kk) {
    sum.all[i,] <- aa[[i]][1,]
    sum.matched[i,] <- aa[[i]][2,]
  }
  if (interactions) {
    sum.all.int <- sum.matched.int <- matrix(NA_real_, nrow = kk*(kk+1)/2, ncol = ncol(aa[[1]]), dimnames = list(NULL, colnames(aa[[1]])))
    to.remove <- rep(FALSE, nrow(sum.all.int))
    int.names <- character(nrow(sum.all.int))
    k <- 1
    for (i in 1:kk) {
      for (j in i:kk) {
        x2 <- X[,i] * X[,j]
        if (all(abs(x2) < sqrt(.Machine$double.eps)) ||
            all(abs(x2 - X[,i]) < sqrt(.Machine$double.eps))) { #prevent interactions within same factors
          to.remove[k] <- TRUE
        }
        else {
          jqoi <- qoi(x2, tt = treat, ww = weights, standardize = standardize)
          sum.all.int[k,] <- jqoi[1,]
          sum.matched.int[k,] <- jqoi[2,]
          if (i == j) {
            int.names[k] <- paste0(nam[i], "\u00B2")
          }
          else {
            int.names[k] <- paste(nam[i], nam[j], sep=" * ")
          }
        }
        k <- k + 1
      }
    }
    rownames(sum.all.int) <- rownames(sum.matched.int) <- int.names

    sum.all <- rbind(sum.all, sum.all.int[!to.remove,,drop = FALSE])
    sum.matched <- rbind(sum.matched, sum.matched.int[!to.remove,,drop = FALSE])
  }

  if (!is.null(object$distance)) {
    ad <- qoi(object$distance, tt = treat, ww = weights, standardize = standardize)
    sum.all <- rbind(ad[1,], sum.all)
    sum.matched <- rbind(ad[2,], sum.matched)
    rownames(sum.all)[1] <- rownames(sum.matched)[1] <- "distance"
  }

  ## Imbalance Reduction
  stat.all <- abs(sum.all[,-(1:2)])
  stat.matched <- abs(sum.matched[,-(1:2)])
  reduction <- 100*(stat.all-stat.matched)/stat.all

  reduction[stat.all == 0 & stat.matched == 0] <- 0
  reduction[stat.all == 0 & stat.matched > 0] <- -Inf

  ## Sample size
  qn <- table(treat[!object$discarded], subclass[!object$discarded])
  dimnames(qn) <- list(c("Control", "Treated"), subclasses)

  small.subclass.control <- subclasses[qn["Control", as.character(subclasses)] <= 1]
  if (length(small.subclass.control) > 0) {
    if (length(small.subclass.control) == 1) warning(paste0("Not enough control units in subclass ", small.subclass.control, "."), call.= FALSE)
    else warning(paste0("Not enough control units in subclasses ", word_list(small.subclass.control), "."), call.= FALSE)
  }

  small.subclass.treated <- subclasses[qn["Treated", as.character(subclasses)] <= 1]
  if (length(small.subclass.treated) > 0) {
    if (length(small.subclass.treated) == 1) warning(paste0("Not enough treated units in subclass ", small.subclass.treated, "."), call.= FALSE)
    else warning(paste0("Not enough treated units in subclasses ", word_list(small.subclass.treated), "."), call.= FALSE)
  }

  if (any(object$discarded)) {
    qn <- cbind(qn, table(treat[object$discarded]))
    colnames(qn)[ncol(qn)] <- "Discarded"
  }
  qn <- rbind(qn, colSums(qn))
  rownames(qn)[nrow(qn)] <- "Total"

  qn <- cbind(qn, rowSums(qn))
  colnames(qn)[ncol(qn)] <- "All"

  ## output
  res <- list(call=object$call, sum.all = sum.all, sum.across = sum.matched,
              sum.subclass = sum.subclass, reduction = reduction,
              qn = qn)
  class(res) <- c("summary.matchit.subclass", "summary.matchit")
  return(res)
}