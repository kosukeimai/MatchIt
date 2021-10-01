### PLOT METHODS-------------------------------------------

plot.matchit <- function(x, type = "qq", interactive = TRUE, which.xs = NULL, ...) {

  type <- tolower(type)
  type <- match_arg(type, c("qq", "ecdf", "density", "jitter", "histogram"))

  if (type %in% c("qq", "ecdf", "density")) {
    matchit.covplot(x, type = type, interactive=interactive,
                    which.xs = which.xs, ...)
  }
  else if (type == "jitter") {
    if (is.null(x$distance)) {
      stop("type = \"jitter\" cannot be used if a distance measure is not estimated or supplied. No plots generated.", call. = FALSE)
    }
    jitter.pscore(x, interactive = interactive,...)
  }
  else if (type =="histogram") {
    if (is.null(x$distance)) {
      stop("type = \"hist\" cannot be used if a distance measure is not estimated or supplied. No plots generated.", call. = FALSE)
    }
    hist.pscore(x,...)
  }
}

plot.matchit.subclass <- function(x, type = "qq", interactive = TRUE, which.xs = NULL, subclass, ...) {
  choice.menu <- function(choices, question) {
    k <- length(choices)-1
    Choices <- data.frame(choices)
    row.names(Choices) <- 0:k
    names(Choices) <- "Choices"
    print.data.frame(Choices, right=FALSE)
    ans <- readline(question)
    while (!ans %in% 0:k) {
      message("Not valid -- please pick one of the choices")
      print.data.frame(Choices, right=FALSE)
      ans <- readline(question)
    }
    return(ans)
  }

  type <- tolower(type)
  type <- match_arg(type, c("qq", "ecdf", "density", "jitter", "histogram"))

  if (type %in% c("qq", "ecdf", "density")) {
    #If subclass = T, index, or range, display all or range of subclasses, using interactive to advance
    #If subclass = F, display aggregate across subclass, using interactive to advance
    #If subclass = NULL, if interactive, use to choose subclass, else display aggregate across subclass

    subclasses <- levels(x$subclass)
    miss.sub <- missing(subclass) || is.null(subclass)
    if (miss.sub || isFALSE(subclass)) which.subclass <- NULL
    else if (isTRUE(subclass)) which.subclass <- subclasses
    else if (!is.atomic(subclass) || !all(subclass %in% seq_along(subclasses))) {
      stop("'subclass' should be TRUE, FALSE, or a vector of subclass indices for which subclass balance is to be displayed.",
           call. = FALSE)
    }
    else which.subclass <- subclasses[subclass]

    if (!is.null(which.subclass)) {
      matchit.covplot.subclass(x, type = type, which.subclass = which.subclass,
                               interactive = interactive, which.xs = which.xs, ...)
    }
    else if (interactive && miss.sub) {
      subclasses <- levels(x$subclass)
      choices <- c("No (Exit)", paste0("Yes: Subclass ", subclasses), "Yes: In aggregate")
      plot.name <- switch(type, "qq" = "quantile-quantile", "ecdf" = "empirical CDF", "density" = "density")
      question <- paste("Would you like to see", plot.name, "plots of any subclasses? ")
      ans <- -1
      while(ans != 0) {
        ans <- as.numeric(choice.menu(choices, question))
        if (ans %in% seq_along(subclasses) && any(x$subclass == subclasses[ans])) {
          matchit.covplot.subclass(x, type = type, which.subclass = subclasses[ans],
                                   interactive = interactive, which.xs = which.xs, ...)
        }
        else if (ans != 0) {
          matchit.covplot(x, type = type, interactive = interactive, which.xs = which.xs, ...)
        }
      }
    }
    else {
      matchit.covplot(x, type = type, interactive = interactive, which.xs = which.xs, ...)
    }
  }
  else if (type=="jitter") {
    if (is.null(x$distance)) {
      stop("type = \"jitter\" cannot be used when no distance variable was estimated or supplied.", call. = FALSE)
    }
    jitter.pscore(x, interactive = interactive, ...)
  }
  else if (type == "histogram") {
    if (is.null(x$distance)) {
      stop("type = \"histogram\" cannot be used when no distance variable was estimated or supplied.", call. = FALSE)
    }
    hist.pscore(x,...)
  }
  invisible(x)
}

plot.summary.matchit <- function(x, abs = TRUE, var.order = "data", threshold = c(.1, .05), position = "bottomright", ...) {

  .pardefault <- par(no.readonly = TRUE)
  on.exit(par(.pardefault))

  sub <- inherits(x, "summary.matchit.subclass")
  matched <- sub || !is.null(x[["sum.matched"]])
  un <- !is.null(x[["sum.all"]])

  standard.sum <- if (un) x[["sum.all"]] else x[[if (sub) "sum.across" else "sum.matched"]]

  if (!"Std. Mean Diff." %in% colnames(standard.sum)) {
    stop("Not appropriate for unstandardized summary.  Run summary() with the standardize = TRUE option, and then plot.", call. = FALSE)
  }

  if (un) {
    sd.all <- x[["sum.all"]][,"Std. Mean Diff."]
  }
  if (matched) {
    sd.matched <- x[[if (sub) "sum.across" else "sum.matched"]][,"Std. Mean Diff."]
  }

  var.names <- rownames(standard.sum)

  var.order <- match_arg(var.order, c("data", "matched", "unmatched", "alphabetical"))

  if (!un && var.order == "unmatched") stop("'var.order' cannot be \"unmatched\" if un = TRUE in the call to summary().", call. = FALSE)
  if (!matched && var.order == "matched") stop("'var.order' cannot be \"matched\" if method = NULL in the original call to matchit().", call. = FALSE)

  if (abs) {
    if (un) sd.all <- abs(sd.all)
    if (matched) sd.matched <- abs(sd.matched)
    xlab <- "Absolute Standardized\nMean Difference"
  }
  else {
    xlab <- "Standardized Mean Difference"
  }

  ord <- switch(var.order,
                "data" = rev(seq_along(var.names)),
                "matched" = order(sd.matched),
                "unmatched" = order(sd.all),
                "alphabetical" = order(var.names, decreasing = TRUE))

  dotchart(if (un) sd.all[ord] else sd.matched[ord],
           labels = var.names[ord], xlab = xlab,
           bg = NA, color = NA, ...)
  abline(v = 0)

  if (sub && length(x$sum.subclass) > 0) {
    for (i in seq_along(x$sum.subclass)) {
      sd.sub <- x$sum.subclass[[i]][,"Std. Mean Diff."]
      if (abs) sd.sub <- abs(sd.sub)
      points(x = sd.sub[ord], y = seq_along(sd.sub),
             pch = as.character(i), col = "gray60", cex = .6)
    }
  }

  if (un) {
    points(x = sd.all[ord], y = seq_along(sd.all),
           pch = 21, bg = "white", col = "black")
  }
  if (matched) {
    points(x = sd.matched[ord], y = seq_along(sd.matched),
           pch = 21, bg = "black", col = "black")
  }

  if (!is.null(threshold)) {
    if (abs) {
      abline(v = threshold, lty = seq_along(threshold))
    }
    else {
      abline(v = threshold, lty = seq_along(threshold))
      abline(v = -threshold, lty = seq_along(threshold))
    }
  }

  if (sum(matched, un) > 1 && !is.null(position)) {
    position <- match_arg(position, c("bottomright", "bottom", "bottomleft", "left",
                                      "topleft", "top", "topright", "right", "center"))
    legend(position, legend = c("All", "Matched"),
           pt.bg = c("white", "black"), pch = 21,
           inset = .015, xpd = TRUE)
  }
  invisible(x)
}

### PRINT METHODS------------------------------------------

print.matchit <- function(x, ...) {
  info <- x[["info"]]
  cal <- !is.null(x[["caliper"]])
  dis <- c("both", "control", "treat")[pmatch(info$discard, c("both", "control", "treat"), 0L)]
  disl <- length(dis) > 0
  nm <- is.null(x[["method"]])
  cat("A matchit object")
  cat(paste0("\n - method: ", info.to.method(info)))
  if (!is.null(x[["distance"]]) || info$mahalanobis || identical(info$distance, "user")) {
    cat("\n - distance: ")
    if (info$mahalanobis) cat("Mahalanobis")
    if (info$distance != "mahalanobis") {
      if (info$mahalanobis) cat(" [matching]\n             ")

      if (info$distance != "user") cat("Propensity score")
      else if (info$distance_is_matrix) cat("User-defined (matrix)")
      else cat("User-defined")

      if (cal || disl) {
        cal.ps <- "" %in% names(x[["caliper"]])
        cat(" [")
        cat(paste(c("matching", "subclassification", "caliper", "common support")[c(!nm && !info$mahalanobis && info$method != "subclass", !nm && info$method == "subclass", cal.ps, disl)], collapse = ", "))
        cat("]")
      }
      if (info$distance != "user") {
        cat("\n             - estimated with ")
        cat(info.to.distance(info))
        if (!is.null(x[["s.weights"]])) {
          if (isTRUE(attr(x[["s.weights"]], "in_ps")))
            cat("\n             - sampling weights included in estimation")
          else cat("\n             - sampling weights not included in estimation")
        }
      }
    }
  }
  if (cal) {
    cat(paste0("\n - caliper: ", paste(vapply(seq_along(x[["caliper"]]), function(z) paste0(if (names(x[["caliper"]])[z] == "") "<distance>" else names(x[["caliper"]])[z],
                                                                                            " (", format(round(x[["caliper"]][z], 3)), ")"), character(1L)),
                                       collapse = ", ")))
  }
  if (disl) {
    cat("\n - common support: ")
    if (dis == "both") cat("units from both groups")
    else if (dis == "treat") cat("treated units")
    else if (dis == "control") cat("control units")
    cat(" dropped")
  }
  cat(paste0("\n - number of obs.: ", length(x[["treat"]]), " (original)", if (!all(x[["weights"]] == 1)) paste0(", ", sum(x[["weights"]] != 0), " (matched)")))
  if (!is.null(x[["s.weights"]])) cat("\n - sampling weights: present")
  if (!is.null(x[["estimand"]])) cat(paste0("\n - target estimand: ", x[["estimand"]]))
  if (!is.null(x[["X"]])) cat(paste0("\n - covariates: ", ifelse(length(names(x[["X"]])) > 40, "too many to name", paste(names(x[["X"]]), collapse = ", "))))
  cat("\n")
  invisible(x)
}

print.summary.matchit <- function(x, digits = max(3, getOption("digits") - 3), ...){

  if (!is.null(x$call)) cat("\nCall:", deparse(x$call), sep = "\n")

  if (!is.null(x$sum.all)) {
    cat("\nSummary of Balance for All Data:\n")
    print.data.frame(round_df_char(x$sum.all[,-7, drop = FALSE], digits, pad = "0", na_vals = "."))
    cat("\n")
  }

  if (!is.null(x$sum.matched)) {
    cat("\nSummary of Balance for Matched Data:\n")
    if (all(is.na(x$sum.matched[,7]))) x$sum.matched <- x$sum.matched[,-7,drop = FALSE] #Remove pair dist if empty
    print.data.frame(round_df_char(x$sum.matched, digits, pad = "0", na_vals = "."))
  }
  if (!is.null(x$reduction)) {
    cat("\nPercent Balance Improvement:\n")
    print.data.frame(round_df_char(x$reduction[,-5, drop = FALSE], 1, pad = "0", na_vals = "."))
  }
  if (!is.null(x$nn)) {
    cat("\nSample Sizes:\n")
    nn <- x$nn
    if (isTRUE(all.equal(nn["All (ESS)",], nn["All",]))) {
      #Don't print ESS if same as full SS
      nn <- nn[rownames(nn) != "All (ESS)",,drop = FALSE]
    }
    if (isTRUE(all.equal(nn["Matched (ESS)",], nn["Matched",]))) {
      #Don't print ESS if same as matched SS
      nn <- nn[rownames(nn) != "Matched (ESS)",,drop = FALSE]
    }
    print.data.frame(round_df_char(nn, 2, pad = " ", na_vals = "."))
  }
  cat("\n")
  invisible(x)
}

print.summary.matchit.subclass <- function(x, digits = max(3, getOption("digits") -  3), ...){

  if (!is.null(x$call)) cat("\nCall:", deparse(x$call), sep = "\n")

  if (!is.null(x$sum.all)) {
    cat("\nSummary of Balance for All Data:\n")
    print.data.frame(round_df_char(x$sum.all[,-7, drop = FALSE], digits, pad = "0", na_vals = "."))
  }

  if (length(x$sum.subclass) > 0) {
    cat("\nSummary of Balance by Subclass:\n")
    for (s in seq_along(x$sum.subclass)) {
      cat(paste0("\n- ", names(x$sum.subclass)[s], "\n"))
      print.data.frame(round_df_char(x$sum.subclass[[s]][,-7, drop = FALSE], digits, pad = "0", na_vals = "."))
    }
    if (!is.null(x$qn)) {
      cat("\nSample Sizes by Subclass:\n")
      print.data.frame(round_df_char(x$qn, 2, pad = " ", na_vals = "."))
    }
  }
  else {
    if (!is.null(x$sum.across)) {
      cat("\nSummary of Balance Across Subclasses\n")
      if (all(is.na(x$sum.across[,7]))) x$sum.across <- x$sum.across[,-7,drop = FALSE]
      print.data.frame(round_df_char(x$sum.across, digits, pad = "0", na_vals = "."))
    }
    if (!is.null(x$reduction)) {
      cat("\nPercent Balance Improvement:\n")
      print.data.frame(round_df_char(x$reduction[,-5, drop = FALSE], 1, pad = "0", na_vals = "."))
    }

    if (!is.null(x$nn)) {
      cat("\nSample Sizes:\n")
      nn <- x$nn
      if (isTRUE(all.equal(nn["All (ESS)",], nn["All",]))) {
        #Don't print ESS if same as full SS
        nn <- nn[rownames(nn) != "All (ESS)",,drop = FALSE]
      }
      if (isTRUE(all.equal(nn["Matched (ESS)",], nn["Matched",]))) {
        #Don't print ESS if same as matched SS
        nn <- nn[rownames(nn) != "Matched (ESS)",,drop = FALSE]
      }
      print.data.frame(round_df_char(nn, 2, pad = " ", na_vals = "."))
    }
  }
  cat("\n")
}

### SUMMARY METHODS----------------------------------------

summary.matchit <- function(object, interactions = FALSE,
                            addlvariables = NULL, standardize = TRUE,
                            data = NULL, pair.dist = TRUE, un = TRUE, improvement = TRUE, ...) {

  #Create covariate matrix; include caliper, exact, and mahvars

  if (is.null(object$X)) {
    X <- matrix(nrow = length(object$treat), ncol = 0)
  }
  else {
    X <- get.covs.matrix(data = object$X)
  }
  X_assign <- get_assign(X)

  if (!is.null(addlvariables)) {
    if (is.character(addlvariables)) {
      if (!is.null(data) && is.data.frame(data)) {
        if (all(addlvariables %in% names(data))) {
          addlvariables <- data[addlvariables]
        }
        else {
          stop("All variables in 'addlvariables' must be in 'data'.", call. = FALSE)
        }
      }
      else {
        stop("If 'addlvariables' is specified as a string, a data frame argument must be supplied to 'data'.", call. = FALSE)
      }
    }
    else if (inherits(addlvariables, "formula")) {
      vars.in.formula <- all.vars(addlvariables)
      if (!is.null(data) && is.data.frame(data)) data <- data.frame(data[names(data) %in% vars.in.formula], object$X[names(data) %in% setdiff(vars.in.formula, names(data))])
      else data <- object$X

      addlvariables <- get.covs.matrix(addlvariables, data = data)
    }
    else if (!is.matrix(addlvariables) && !is.data.frame(addlvariables)) {
      stop("The argument to 'addlvariables' must be in one of the accepted forms. See ?summary.matchit for details.", call. = FALSE)
    }

    if (is.data.frame(addlvariables)) {
      addlvariables <- get.covs.matrix(data = addlvariables)
    }

    addl_assign <- get_assign(addlvariables)
    X <- cbind(X, addlvariables[, setdiff(colnames(addlvariables), colnames(X)), drop = FALSE])
  }

  treat <- object$treat
  weights <- object$weights
  s.weights <- if (is.null(object$s.weights)) rep(1, length(weights)) else object$s.weights

  kk <- ncol(X)

  if (kk > 0) {
    nam <- colnames(X)
    nam[startsWith(nam, "`") & endsWith(nam, "`")] <- substr(nam[startsWith(nam, "`") & endsWith(nam, "`")],
                                                             2, nchar(nam[startsWith(nam, "`") & endsWith(nam, "`")]) - 1)
  }

  matched <- !is.null(object$info$method)
  un <- un || !matched

  if (standardize) {
    s.d.denom <- switch(object$estimand,
                        "ATT" = "treated",
                        "ATC" = "control",
                        "ATE" = "pooled")
  }
  else s.d.denom <- NULL

  ## Summary Stats
  if (kk > 0) {
    if (un) {
      aa.all <- setNames(lapply(seq_len(kk), function(i) qoi(X[,i], tt = treat, ww = NULL, s.weights = s.weights,
                                                             standardize = standardize, s.d.denom = s.d.denom)),
                         colnames(X))
      sum.all <- do.call("rbind", aa.all)
      dimnames(sum.all) <- list(nam, names(aa.all[[1]]))

      sum.all.int <- NULL
    }

    if (matched) {
      aa.matched <- setNames(lapply(seq_len(kk), function(i) qoi(X[,i], tt = treat, ww = weights, s.weights = s.weights,
                                                                 subclass = object$subclass, mm = object$match.matrix,
                                                                 standardize = standardize, s.d.denom = s.d.denom,
                                                                 compute.pair.dist = pair.dist)),
                             colnames(X))
      sum.matched <- do.call("rbind", aa.matched)
      dimnames(sum.matched) <- list(nam, names(aa.matched[[1]]))

      sum.matched.int <- NULL
    }

    if (interactions) {
      n.int <- kk*(kk+1)/2
      if (un) sum.all.int <- matrix(NA_real_, nrow = n.int, ncol = length(aa.all[[1]]), dimnames = list(NULL, names(aa.all[[1]])))
      if (matched) sum.matched.int <- matrix(NA_real_, nrow = n.int, ncol = length(aa.matched[[1]]), dimnames = list(NULL, names(aa.matched[[1]])))

      to.remove <- rep(FALSE, n.int)
      int.names <- character(n.int)
      k <- 1
      for (i in 1:kk) {
        for (j in i:kk) {
          x2 <- X[,i] * X[,j]
          if (all(abs(x2) < sqrt(.Machine$double.eps)) ||
              all(abs(x2 - X[,i]) < sqrt(.Machine$double.eps))) { #prevent interactions within same factors
            to.remove[k] <- TRUE
          }
          else {
            if (un) {
              sum.all.int[k,] <- qoi(x2, tt = treat, ww = NULL, s.weights = s.weights,
                                     standardize = standardize, s.d.denom = s.d.denom)
            }
            if (matched) {
              sum.matched.int[k,] <- qoi(x2, tt = treat, ww = weights, s.weights = s.weights,
                                         subclass = object$subclass, mm = object$match.matrix,
                                         standardize = standardize, s.d.denom = s.d.denom,
                                         compute.pair.dist = pair.dist)
            }
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

      if (un) {
        rownames(sum.all.int) <- int.names
        sum.all <- rbind(sum.all, sum.all.int[!to.remove,,drop = FALSE])
      }
      if (matched) {
        rownames(sum.matched.int) <- int.names
        sum.matched <- rbind(sum.matched, sum.matched.int[!to.remove,,drop = FALSE])
      }
    }
  }

  if (!is.null(object$distance)) {
    if (un) {
      ad.all <- qoi(object$distance, tt = treat, ww = NULL, s.weights = s.weights,
                    standardize = standardize, s.d.denom = s.d.denom)
      if (!exists("sum.all", inherits = FALSE)) {
        sum.all <- matrix(ad.all, nrow = 1, dimnames = list("distance", names(ad.all)))
      }
      else {
        sum.all <- rbind(ad.all, sum.all)
        rownames(sum.all)[1] <- "distance"
      }
    }
    if (matched) {
      ad.matched <- qoi(object$distance, tt = treat, ww = weights, s.weights = s.weights,
                        subclass = object$subclass, mm = object$match.matrix, standardize = standardize,
                        s.d.denom = s.d.denom, compute.pair.dist = pair.dist)
      if (!exists("sum.matched", inherits = FALSE)) {
        sum.matched <- matrix(ad.matched, nrow = 1, dimnames = list("distance", names(ad.matched)))
      }
      else {
        sum.matched <- rbind(ad.matched, sum.matched)
        rownames(sum.matched)[1] <- "distance"
      }
    }
  }

  ## Imbalance Reduction
  if (matched && un && improvement) {
    reduction <- matrix(NA_real_, nrow = nrow(sum.all), ncol = ncol(sum.all) - 2,
                        dimnames = list(rownames(sum.all), colnames(sum.all)[-(1:2)]))
    stat.all <- abs(sum.all[,-(1:2), drop = FALSE])
    stat.matched <- abs(sum.matched[,-(1:2), drop = FALSE])

    #Everything but variance ratios
    reduction[,-2] <- 100*(stat.all[,-2]-stat.matched[,-2])/stat.all[,-2]

    #Just variance ratios; turn to log first
    vr.all <- abs(log(stat.all[,2]))
    vr.matched <- abs(log(stat.matched[,2]))
    reduction[,2] <- 100*(vr.all-vr.matched)/vr.all

    reduction[stat.all == 0 & stat.matched == 0] <- 0
    reduction[stat.all == 0 & stat.matched > 0] <- -Inf
  }
  else {
    reduction <- NULL
  }

  ## output
  res <- list(call = object$call, nn = object$nn, sum.all = if (un) sum.all,
              sum.matched = if (matched) sum.matched, reduction = reduction)
  class(res) <- "summary.matchit"
  return(res)
}

summary.matchit.subclass <- function(object, interactions = FALSE,
                                     addlvariables = NULL, standardize = TRUE,
                                     data = NULL, pair.dist = FALSE,
                                     subclass = FALSE, un = TRUE, improvement = TRUE, ...) {

  #Create covariate matrix
  X <- get.covs.matrix(data = object$X)

  if (!is.null(addlvariables)) {
    if (is.character(addlvariables)) {
      if (!is.null(data) && is.data.frame(data)) {
        if (all(addlvariables %in% names(data))) {
          addlvariables <- data[addlvariables]
        }
        else {
          stop("All variables in 'addlvariables' must be in 'data'.", call. = FALSE)
        }
      }
      else {
        stop("If 'addlvariables' is specified as a string, a data frame argument must be supplied to 'data'.", call. = FALSE)
      }
    }
    else if (inherits(addlvariables, "formula")) {
      vars.in.formula <- all.vars(addlvariables)
      if (!is.null(data) && is.data.frame(data)) data <- data.frame(data[names(data) %in% vars.in.formula],
                                                                    object$X[names(data) %in% setdiff(vars.in.formula, names(data))])
      else data <- object$X

      addlvariables <- get.covs.matrix(addlvariables, data = data)
    }
    else if (!is.matrix(addlvariables) && !is.data.frame(addlvariables)) {
      stop("The argument to 'addlvariables' must be in one of the accepted forms. See ?summary.matchit for details.", call. = FALSE)
    }

    if (is.data.frame(addlvariables)) {
      if (!all(vapply(addlvariables, is.numeric, logical(1L)))) {
        addlvariables <- get.covs.matrix(data = addlvariables)
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
  s.weights <- if (is.null(object$s.weights)) rep(1, length(weights)) else object$s.weights
  subclass <- object$subclass

  nam <- colnames(X)

  kk <- ncol(X)
  subclasses <- levels(subclass)

  if (standardize) {
    s.d.denom <- switch(object$estimand,
                        "ATT" = "treated",
                        "ATC" = "control",
                        "ATE" = "pooled")
  }
  else s.d.denom <- NULL

  if (isTRUE(which.subclass)) which.subclass <- subclasses
  else if (isFALSE(which.subclass)) which.subclass <- NULL
  else if (!is.atomic(which.subclass) || !all(which.subclass %in% seq_along(subclasses))) {
    stop("'subclass' should be TRUE, FALSE, or a vector of subclass indices for which subclass balance is to be displayed.")
  }
  else which.subclass <- subclasses[which.subclass]

  matched <- TRUE #always compute aggregate balance so plot.summary can use it
  subs <- !is.null(which.subclass)

  ## Aggregate Subclass
  #Use the estimated weights to compute aggregate balance.
  ## Summary Stats

  sum.all <- sum.matched <- sum.subclass <- reduction <- NULL

  if (un) {
    aa.all <- setNames(lapply(seq_len(kk), function(i) qoi(X[,i], tt = treat, ww = NULL, s.weights = s.weights,
                                                           standardize = standardize, s.d.denom = s.d.denom)),
                       colnames(X))
    sum.all <- do.call("rbind", aa.all)
    dimnames(sum.all) <- list(nam, names(aa.all[[1]]))

    sum.all.int <- NULL
  }

  if (matched) {
    aa.matched <- setNames(lapply(seq_len(kk), function(i) qoi(X[,i], tt = treat, ww = weights, s.weights = s.weights,
                                                               subclass = subclass, standardize = standardize,
                                                               s.d.denom = s.d.denom, compute.pair.dist = pair.dist)),
                           colnames(X))
    sum.matched <- do.call("rbind", aa.matched)
    dimnames(sum.matched) <- list(nam, names(aa.matched[[1]]))

    sum.matched.int <- NULL
  }

  if (interactions) {
    n.int <- kk*(kk+1)/2
    if (un) sum.all.int <- matrix(NA_real_, nrow = n.int, ncol = length(aa.all[[1]]), dimnames = list(NULL, names(aa.all[[1]])))
    if (matched) sum.matched.int <- matrix(NA_real_, nrow = n.int, ncol = length(aa.matched[[1]]), dimnames = list(NULL, names(aa.matched[[1]])))

    to.remove <- rep(FALSE, n.int)
    int.names <- character(n.int)
    k <- 1
    for (i in 1:kk) {
      for (j in i:kk) {
        x2 <- X[,i] * X[,j]
        if (all(abs(x2) < sqrt(.Machine$double.eps)) ||
            all(abs(x2 - X[,i]) < sqrt(.Machine$double.eps))) { #prevent interactions within same factors
          to.remove[k] <- TRUE
        }
        else {
          if (un) {
            sum.all.int[k,] <- qoi(x2, tt = treat, ww = NULL, s.weights = s.weights,
                                   standardize = standardize, s.d.denom = s.d.denom)
          }
          if (matched) {
            sum.matched.int[k,] <- qoi(x2, tt = treat, ww = weights, s.weights = s.weights,
                                       subclass = subclass, standardize = standardize,
                                       compute.pair.dist = pair.dist)
          }
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

    if (un) {
      rownames(sum.all.int) <- int.names
      sum.all <- rbind(sum.all, sum.all.int[!to.remove,,drop = FALSE])
    }
    if (matched) {
      rownames(sum.matched.int) <- int.names
      sum.matched <- rbind(sum.matched, sum.matched.int[!to.remove,,drop = FALSE])
    }
  }

  if (!is.null(object$distance)) {
    if (un) {
      ad.all <- qoi(object$distance, tt = treat, ww = NULL, s.weights = s.weights,
                    standardize = standardize, s.d.denom = s.d.denom)
      sum.all <- rbind(ad.all, sum.all)
      rownames(sum.all)[1] <- "distance"
    }
    if (matched) {
      ad.matched <- qoi(object$distance, tt = treat, ww = weights, s.weights = s.weights,
                        subclass = subclass, standardize = standardize,
                        s.d.denom = s.d.denom, compute.pair.dist = pair.dist)
      sum.matched <- rbind(ad.matched, sum.matched)
      rownames(sum.matched)[1] <- "distance"
    }
  }

  ## Imbalance Reduction
  if (un && matched && improvement) {
    stat.all <- abs(sum.all[,-(1:2)])
    stat.matched <- abs(sum.matched[,-(1:2)])
    reduction <- 100*(stat.all-stat.matched)/stat.all

    reduction[stat.all == 0 & stat.matched == 0] <- 0
    reduction[stat.all == 0 & stat.matched > 0] <- -Inf
  }

  ## By Subclass
  if (subs) {
    sum.subclass <- lapply(which.subclass, function(s) {

      #qoi.subclass only returns unmatched stats, which is all we need within
      #subclasses. Otherwise, identical to matched stats.
      aa <- setNames(lapply(seq_len(kk), function(i) {
        qoi.subclass(X[,i], tt = treat, s.weights = s.weights, subclass = subclass, s.d.denom = s.d.denom, standardize = standardize, which.subclass = s)
      }), colnames(X))

      sum.sub <- matrix(NA_real_, nrow = kk, ncol = ncol(aa[[1]]), dimnames = list(nam, colnames(aa[[1]])))

      sum.sub.int <- NULL
      for (i in 1:kk) {
        sum.sub[i,] <- aa[[i]]
      }
      if (interactions) {
        sum.sub.int <- matrix(NA_real_, nrow = kk*(kk+1)/2, ncol = length(aa[[1]]), dimnames = list(NULL, names(aa[[1]])))
        to.remove <- rep(FALSE, nrow(sum.sub.int))
        int.names <- character(nrow(sum.sub.int))
        k <- 1
        for (i in 1:kk) {
          for (j in i:kk) {
            if (!to.remove[k]) { #to.remove defined above
              x2 <- X[,i] * X[,j]
              jqoi <- qoi.subclass(x2, tt = treat, s.weights = s.weights, subclass = subclass, s.d.denom = s.d.denom, standardize = standardize, which.subclass = s)
              sum.sub.int[k,] <- jqoi
              if (i == j) {
                int.names[k] <- paste0(nam[i], "\u00B2")
              }
              else {
                int.names[k] <- paste(nam[i], nam[j], sep = " * ")
              }
            }
            k <- k + 1
          }
        }
        rownames(sum.sub.int) <- int.names

        sum.sub <- rbind(sum.sub, sum.sub.int[!to.remove,,drop = FALSE])
      }

      if (!is.null(object$distance)) {
        ad <- qoi.subclass(object$distance, tt = treat, s.weights = s.weights, subclass = subclass,
                           s.d.denom = s.d.denom, standardize = standardize, which.subclass = s)
        sum.sub <- rbind(ad, sum.sub)
        rownames(sum.sub)[1] <- "distance"
      }

      return(sum.sub)
    })
    names(sum.subclass) <- paste("Subclass", which.subclass)
  }

  ## Sample size
  qn <- object$qn

  if (subs) {
    small.subclass.control <- which.subclass[qn["Control", as.character(which.subclass)] <= 1]
    if (length(small.subclass.control) > 0) {
      if (length(small.subclass.control) == 1) warning(paste0("Not enough control units in subclass ", small.subclass.control, "."), call.= FALSE)
      else warning(paste0("Not enough control units in subclasses ", word_list(small.subclass.control), "."), call.= FALSE)
    }

    small.subclass.treated <- which.subclass[qn["Treated", as.character(which.subclass)] <= 1]
    if (length(small.subclass.treated) > 0) {
      if (length(small.subclass.treated) == 1) warning(paste0("Not enough treated units in subclass ", small.subclass.treated, "."), call.= FALSE)
      else warning(paste0("Not enough treated units in subclasses ", word_list(small.subclass.treated), "."), call.= FALSE)
    }
  }

  ## output
  res <- list(call=object$call, sum.all = sum.all, sum.across = sum.matched,
              sum.subclass = sum.subclass, reduction = reduction,
              qn = qn, nn = object$nn)
  class(res) <- c("summary.matchit.subclass", "summary.matchit")
  return(res)
}