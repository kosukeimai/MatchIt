### PLOT METHODS-------------------------------------------

plot.matchit <- function(x, type = "qq", interactive = TRUE, which.xs = NULL, ...) {

  type <- tolower(type)
  type <- match_arg(type, c("qq", "jitter", "histogram", "ecdf"))

  if (type == "qq") {
    matchit.qqplot(x, interactive=interactive,
                   which.xs = which.xs, ...)
  }
  else if (type == "ecdf") {
    matchit.ecdfplot(x, interactive=interactive,
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

plot.matchit.subclass <- function(x, type = "qq", interactive = TRUE, which.xs = NULL, subclass = NULL, ...) {
  choice.menu <- function(choices, question) {
    k <- length(choices)-1
    Choices <- data.frame(choices)
    row.names(Choices) <- 0:k
    names(Choices) <- "Choices"
    print.data.frame(Choices, right=FALSE)
    ans <- readline(question)
    while(!ans %in% 0:k) {
      print("Not valid -- please pick one of the choices")
      print.data.frame(Choices, right=FALSE)
      ans <- readline(question)
    }
    return(ans)
  }

  type <- tolower(type)
  type <- match_arg(type, c("qq", "ecdf", "jitter", "histogram"))

  if (type == "qq"){
    if (interactive) {
      subclasses <- levels(x$subclass)
      choices <- c("No (Exit)", paste0("Yes: Subclass ", subclasses), "Yes: In aggregate")
      question <- "Would you like to see quantile-quantile plots of any subclasses? "
      ans <- -1
      while(ans != 0) {
        ans <- as.numeric(choice.menu(choices, question))
        if (ans %in% seq_along(subclasses) && sum(x$subclass == subclasses[ans]) > 0) {
          matchit.qqplot(x, which.subclass = subclasses[ans],
                         interactive = interactive, which.xs = which.xs,...)
        }
        else if (ans != 0) {
          matchit.qqplot(x, interactive = interactive, which.xs = which.xs,...)
        }
      }
    }
    else {
      matchit.qqplot(x,  which.subclass = subclass,
                     interactive = interactive, which.xs = which.xs,...)
    }
  }
  else if (type == "ecdf"){
    if (interactive) {
      subclasses <- levels(x$subclass)
      choices <- c("No (Exit)", paste0("Yes: Subclass ", subclasses), "Yes: In aggregate")
      question <- "Would you like to see emprical CDF plots of any subclasses? "
      ans <- -1
      while(ans != 0) {
        ans <- as.numeric(choice.menu(choices, question))
        if (ans %in% seq_along(subclasses) && sum(x$subclass == subclasses[ans]) > 0) {
          matchit.ecdfplot(x, which.subclass = subclasses[ans],
                         interactive = interactive, which.xs = which.xs,...)
        }
        else if (ans != 0) {
          matchit.ecdfplot(x, interactive = interactive, which.xs = which.xs,...)
        }
      }
    }
    else {
      matchit.ecdfplot(x,  which.subclass = subclass,
                       interactive = interactive, which.xs = which.xs,...)
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
    xlab <- "Absolute Standardized\nMean Difference"
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

### PRINT METHODS------------------------------------------

print.matchit <- function(x, ...) {
  info <- x[["info"]]
  cal <- !is.null(x[["caliper"]])
  dis <- c("both", "control", "treat")[pmatch(info$discard, c("both", "control", "treat"), 0L)]
  disl <- length(dis) > 0
  cat("A matchit object\n")
  cat(paste0(" - method: ", info.to.method(info), "\n"))
  if (!is.null(x[["distance"]])) {
    cat(" - distance: ")
    if (info$mahalanobis) cat("Mahalanobis")
    if (info$distance != "mahalanobis") {
      if (info$mahalanobis) cat(" [matching]\n             ")
      if (info$distance == "user") cat("User-defined") else cat("Propensity score")

      if (cal || disl) {
        cat(" [")
        cat(paste(c("matching", "subclassification", "caliper", "common support")[c(!info$mahalanobis && info$method != "subclass", info$method == "subclass", cal, disl)], collapse = ", "))
        cat("]")
      }
      if (info$distance != "user") {
        cat("\n              - estimated with ")
        cat(info.to.distance(info))
      }
    }
    cat("\n")
  }
  if (cal) {
    cat(paste0(" - caliper: ", paste(vapply(seq_along(x[["caliper"]]), function(z) paste0(if (names(x[["caliper"]])[z] == "") "<distance>" else names(x[["caliper"]])[z],
                                                                                          " (", format(round(x[["caliper"]][z], 3)), ")"), character(1L)),
                                     collapse = ", "), "\n"))
  }
  if (disl) {
    cat(" - common support: ")
    if (dis == "both") cat("units from both groups")
    else if (dis == "treat") cat("treated units")
    else if (dis == "control") cat("control units")
    cat(" dropped\n")
  }
  cat(paste0(" - number of obs.: ", length(x[["treat"]]), " (original), ", sum(x[["weights"]] != 0), " (matched)\n"))
  cat(paste0(" - target estimand: ", x[["estimand"]], "\n"))
  if (!is.null(x[["X"]])) cat(paste0(" - covariates: ", ifelse(length(names(x[["X"]])) > 40, "too many to name", paste(names(x[["X"]]), collapse = ", ")), "\n"))
  invisible(x)
}

print.summary.matchit <- function(x, digits = max(3, getOption("digits") - 3), ...){

  cat("\nCall:", deparse(x$call), sep = "\n")
  cat("\nSummary of Balance for All Data:\n")
  print.data.frame(round_df_char(x$sum.all[,-7, drop = FALSE], digits, pad = "0", na_vals = "."))
  cat("\n")

  if(!is.null(x$sum.matched)) {
    cat("\nSummary of Balance for Matched Data:\n")
    print.data.frame(round_df_char(x$sum.matched, digits, pad = "0", na_vals = "."))
    cat("\nPercent Balance Improvement:\n")
    print.data.frame(round_df_char(x$reduction[,-5, drop = FALSE], 1, pad = "0", na_vals = "."))
    cat("\nSample Sizes:\n")
    print.data.frame(round_df_char(x$nn, 1, pad = "0", na_vals = "."))
    cat("\n")
  }
  invisible(x)
}

print.summary.matchit.subclass <- function(x, digits = max(3, getOption("digits") -  3), ...){
  cat("\nCall:", deparse(x$call), sep = "\n")
  cat("\nSummary of Balance for All Data:\n")
  print.data.frame(round_df_char(x$sum.all[,-7, drop = FALSE], digits, pad = "0", na_vals = "."))

  if (length(x$sum.subclass) > 0) {
    cat("\nSummary of Balance by Subclass:\n")
    for (s in seq_along(x$sum.subclass)) {
      cat(paste0("\n- ", names(x$sum.subclass)[s], "\n"))
      print.data.frame(round_df_char(x$sum.subclass[[s]][,-7, drop = FALSE], digits, pad = "0", na_vals = "."))
    }
  }
  cat("\nSummary of Balance Across Subclasses\n")
  print.data.frame(round_df_char(x$sum.across, digits, pad = "0", na_vals = "."))
  cat("\nPercent Balance Improvement:\n")
  print.data.frame(round_df_char(x$reduction[,-5, drop = FALSE], 1, pad = "0", na_vals = "."))
  cat("\nSample Sizes by Subclass:\n")
  print.data.frame(round_df_char(x$qn, 1, pad = "0", na_vals = "."))
  cat("\n")
}



### SUMMARY METHODS----------------------------------------

summary.matchit <- function(object, interactions = FALSE,
                            addlvariables = NULL, standardize = TRUE,
                            data = NULL, ...) {

  #Create covariate matrix; include caliper, exact, and mahvars

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
      if (!is.null(data) && is.data.frame(data)) data <- data.frame(data[names(data) %in% vars.in.formula], object$X[names(data) %in% setdiff(vars.in.formula, names(data))])
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

    X <- cbind(X, addlvariables[, setdiff(colnames(addlvariables), colnames(X)), drop = FALSE])
  }

  treat <- object$treat
  weights <- object$weights
  nam <- colnames(X)
  nam[startsWith(nam, "`") & endsWith(nam, "`")] <- substr(nam[startsWith(nam, "`") & endsWith(nam, "`")],
                                                           2, nchar(nam[startsWith(nam, "`") & endsWith(nam, "`")]) - 1)

  kk <- ncol(X)

  if (standardize) {
    s.d.denom <- switch(object$estimand,
                        "ATT" = "treated",
                        "ATC" = "control",
                        "ATE" = "pooled")
  }
  else s.d.denom <- NULL

  ## Summary Stats
  aa <- setNames(lapply(seq_len(kk), function(i) qoi(X[,i], tt = treat, ww = weights, subclass = object$subclass,
                                                     mm = object$match.matrix, standardize = standardize, s.d.denom = s.d.denom)),
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
          jqoi <- qoi(x2, tt = treat, ww = weights, subclass = object$subclass,
                      mm = object$match.matrix, standardize = standardize, s.d.denom = s.d.denom)
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
    ad <- qoi(object$distance, tt = treat, ww = weights, subclass = object$subclass,
              mm = object$match.matrix, standardize = standardize, s.d.denom = s.d.denom)
    sum.all <- rbind(ad[1,], sum.all)
    sum.matched <- rbind(ad[2,], sum.matched)
    rownames(sum.all)[1] <- rownames(sum.matched)[1] <- "distance"
  }

  ## Imbalance Reduction
  reduction <- matrix(NA_real_, nrow = nrow(sum.all), ncol = ncol(sum.all) - 2,
                      dimnames = list(rownames(sum.all), colnames(sum.all)[-(1:2)]))
  stat.all <- abs(sum.all[,-(1:2)])
  stat.matched <- abs(sum.matched[,-(1:2)])

  #Everything but variance ratios
  reduction[,-2] <- 100*(stat.all[,-2]-stat.matched[,-2])/stat.all[,-2]

  #Just variance ratios; turn to log first
  vr.all <- abs(log(stat.all[,2]))
  vr.matched <- abs(log(stat.matched[,2]))
  reduction[,2] <- 100*(vr.all-vr.matched)/vr.all

  reduction[stat.all == 0 & stat.matched == 0] <- 0
  reduction[stat.all == 0 & stat.matched > 0] <- -Inf

  ## output
  res <- list(call = object$call, nn = object$nn, sum.all = sum.all,
              sum.matched = sum.matched, reduction = reduction)
  class(res) <- "summary.matchit"
  return(res)
}

summary.matchit.subclass <- function(object, interactions = FALSE,
                                     addlvariables = NULL, standardize = TRUE,
                                     data = NULL, subclass = FALSE, ...) {

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
  else if (!is.atomic(which.subclass) || !all(which.subclass %in% subclasses)) {
    stop("'subclass' should be TRUE, FALSE, or a vector of subclass indices for which subclass balance is to be displayed.")
  }
  else which.subclass <- subclasses[which.subclass]

  ## Aggregate Subclass
  #Use the estimated weights to compute aggregate balance.
  ## Summary Stats
  aa <- setNames(lapply(seq_len(kk), function(i) qoi(X[,i], tt = treat, ww = weights, subclass = subclass, standardize = standardize)),
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
          jqoi <- qoi(x2, tt = treat, ww = weights, subclass = subclass, standardize = standardize)
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
    ad <- qoi(object$distance, tt = treat, ww = weights, subclass = subclass, standardize = standardize)
    sum.all <- rbind(ad[1,], sum.all)
    sum.matched <- rbind(ad[2,], sum.matched)
    rownames(sum.all)[1] <- rownames(sum.matched)[1] <- "distance"
  }

  ## By Subclass
  sum.subclass <- lapply(which.subclass, function(s) {

    #qoi without weights only returns unmatched stats, which is all we need within
    #subclasses. Otherwise, identical to matched stats.
    aa <- setNames(lapply(seq_len(kk), function(i) {
      qoi.subclass(X[,i], tt = treat, subclass = subclass, s.d.denom = s.d.denom, standardize = standardize, which.subclass = s)
    }), colnames(X))

    sum.sub <- matrix(NA_real_, nrow = kk, ncol = ncol(aa[[1]]), dimnames = list(nam, colnames(aa[[1]])))

    sum.sub.int <- NULL
    for (i in 1:kk) {
      sum.sub[i,] <- aa[[i]]
    }
    if (interactions) {
      sum.sub.int <- matrix(NA_real_, nrow = kk*(kk+1)/2, ncol = ncol(aa[[1]]), dimnames = list(NULL, colnames(aa[[1]])))
      to.remove <- rep(FALSE, nrow(sum.sub.int))
      int.names <- character(nrow(sum.sub.int))
      k <- 1
      for (i in 1:kk) {
        for (j in i:kk) {
          if (!to.remove[k]) { #to.remove defined above
            x2 <- X[,i] * X[,j]
            jqoi <- qoi.subclass(x2, tt = treat, subclass = subclass, s.d.denom = s.d.denom, standardize = standardize, which.subclass = s)
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
      ad <- qoi.subclass(object$distance, tt = treat, subclass = subclass, s.d.denom = s.d.denom, standardize = standardize, which.subclass = s)
      sum.sub <- rbind(ad, sum.sub)
      rownames(sum.sub)[1] <- "distance"
    }

    return(sum.sub)
  })
  if (!is.null(which.subclass)) names(sum.subclass) <- paste("Subclass", which.subclass)

  ## Imbalance Reduction
  stat.all <- abs(sum.all[,-(1:2)])
  stat.matched <- abs(sum.matched[,-(1:2)])
  reduction <- 100*(stat.all-stat.matched)/stat.all

  reduction[stat.all == 0 & stat.matched == 0] <- 0
  reduction[stat.all == 0 & stat.matched > 0] <- -Inf

  ## Sample size
  qn <- table(treat[!object$discarded], subclass[!object$discarded])
  dimnames(qn) <- list(c("Control", "Treated"), subclasses)

  if (!is.null(which.subclass)) {
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