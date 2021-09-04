cem_matchit <- function(treat, X, cutpoints = "sturges", grouping = list(), k2k = FALSE, k2k.method = "mahalanobis", mpower = 2, estimand = "ATT") {
  #In-house implementation of cem. Basically the same except:
  #treat is a vector if treatment status, not the name of a variable
  #X is a data.frame of covariates
  #when cutpoints are given as integer or string, they define the number of bins, not the number of breakpoints. "ss" is no longer allowed.
  #When k2k = TRUE, subclasses are created for each pair, mimicking true matching, not each covariate combo.
  #k2k.method is used instead of method. When k2k.method = NULL, units are matched based on order rather than random. Default is "mahalanobis" (not available in cem).
  #k2k now works with single covariates (previously it was ignored). k2k uses original variables, not coarsened versions

  if (k2k && !is.null(k2k.method)) {
    X.match <- scale(get.covs.matrix(data = X), center = FALSE)
    k2k.method <- tolower(k2k.method)
    k2k.method <- match_arg(k2k.method, c("mahalanobis", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"))
    if (k2k.method == "mahalanobis") mahSigma_inv <- generalized_inverse(cov(X.match))
  }
  is.numeric.cov <- setNames(vapply(X, is.numeric, logical(1L)), names(X))

  #Process cutpoints
  if (!is.list(cutpoints)) {
    cutpoints <- setNames(as.list(rep(cutpoints, sum(is.numeric.cov))), names(X)[is.numeric.cov])
  }

  if (is.null(names(cutpoints))) stop("'cutpoints' must be a named list of binning values with an element for each numeric variable.", call. = FALSE)
  bad.names <- setdiff(names(cutpoints), names(X))
  nb <- length(bad.names)
  if (nb > 0) {
    warning(paste0(ngettext(nb, "The variable ", "The variables "),
                   word_list(bad.names, quotes = 2), " named in 'cutpoints' ",
                   ngettext(nb, "is ", "are "), "not in the variables supplied to matchit() and will be ignored."),
            immediate. = TRUE, call. = FALSE)
    cutpoints[bad.names] <- NULL
  }

  bad.cuts <- setNames(rep(FALSE, length(cutpoints)), names(cutpoints))
  for (i in names(cutpoints)) {
    if (length(cutpoints[[i]]) == 0) {
      cutpoints[[i]] <- "sturges"
    }
    else if (length(cutpoints[[i]]) == 1) {
      if (is.character(cutpoints[[i]])) {

        bad.cuts[i] <- !(startsWith(cutpoints[[i]], "q") && can_str2num(substring(cutpoints[[i]], 2))) &&
          is.na(pmatch(cutpoints[[i]], c("sturges", "fd", "scott")))
      }
      else if (is.numeric(cutpoints[[i]])) {
        if      (!is.finite(cutpoints[[i]]) || cutpoints[[i]] < 0) bad.cuts[i] <- TRUE
        if      (cutpoints[[i]] == 0) is.numeric.cov[i] <- FALSE #Will not be binned
        else if (cutpoints[[i]] == 1) X[[i]] <- NULL #Removing from X, still in X.match
      }
    }
    else {
      bad.cuts[i] <- !is.numeric(cutpoints[[i]])
    }
  }
  if (any(bad.cuts)) {
    stop(paste0("All entries in the list supplied to 'cutpoints' must be one of the following:",
                "\n\t- a string containing the name of an allowable binning method",
                "\n\t- a single number corresponding to the number of bins",
                "\n\t- a numeric vector containing the cut points separating bins",
                "\nIncorrectly specified ", ngettext(sum(bad.cuts), "variable:\n\t", "variables:\n\t"),
                paste(names(cutpoints)[bad.cuts], collapse = ", ")), call. = FALSE)
  }

  #Process grouping
  if (!is.null(grouping) && !is.null(names(grouping))) {
    X[names(grouping)] <- lapply(names(grouping), function(g) {
      x <- X[[g]]
      groups <- grouping[[g]]

      for (i in seq_along(groups)) {
        x[x %in% groups[[i]]] <- groups[[i]][1]
      }
      x
    })
    cutpoints[names(cutpoints) %in% names(grouping)] <- NULL
  }

  #Create bins for numeric variables
  for (i in names(X)[is.numeric.cov]) {
    if (is.null(cutpoints) || !i %in% names(cutpoints)) bins <- "sturges"
    else bins <- cutpoints[[i]]

    if (is.character(bins)) {
      if (startsWith(bins, "q") || can_str2num(substring(bins, 2))) {
        #Quantile bins
        q <- str2num(substring(bins, 2))
        bins <- quantile(X[[i]], probs = seq(1/q, 1 - 1/q, by = 1/q), names = FALSE) #Outer boundaries will be added later
      }
      else {
        bins <- match_arg(tolower(bins), c("sturges", "fd", "scott"))
        bins <- switch(bins,
                       sturges = nclass.Sturges(X[[i]]),
                       fd = nclass.FD(X[[i]]),
                       scott = nclass.scott(X[[i]]))
        #Breaks is now a single number
      }
    }

    if (length(bins) == 1) {
      #cutpoints is number of bins, unlike in cem
      breaks <- seq(min(X[[i]]), max(X[[i]]), length = bins + 1)
      breaks[c(1, bins + 1)] <- c(-Inf, Inf)
    }
    else {
      breaks <- c(-Inf, sort(unique(bins)), Inf)
    }

    X[[i]] <- findInterval(X[[i]], breaks)
  }

  #Exact match
  xx <- exactify(X, names(treat))
  cc <- do.call("intersect", unname(split(xx, treat)))

  if (length(cc) == 0) {
    stop("No units were matched. Try coarsening the variables further or decrease the number of variables to match on.", call. = FALSE)
  }

  subclass <- setNames(match(xx, cc), names(treat))

  extra.sub <- max(subclass, na.rm = TRUE)

  if (k2k) {

    na.sub <- is.na(subclass)

    s <- switch(estimand, "ATC" = 0, 1)

    for (i in which(tabulateC(subclass[!na.sub]) > 2)) {

      in.sub <- which(!na.sub & subclass == i)

      #Compute distance matrix; all 0s if k2k.method = NULL for matched based on data order
      if (is.null(k2k.method)) dist.mat <- matrix(0, nrow = length(in.sub), ncol = length(in.sub),
                                                  dimnames = list(names(treat)[in.sub], names(treat)[in.sub]))
      else if (k2k.method == "mahalanobis") {
        dist.mat <- matrix(0, nrow = length(in.sub), ncol = length(in.sub),
                           dimnames = list(names(treat)[in.sub], names(treat)[in.sub]))
        for (t in seq_along(in.sub)) dist.mat[t,] <- mahalanobis(X.match[in.sub,,drop = FALSE], X.match[in.sub[t],], mahSigma_inv, TRUE)
      }
      else dist.mat <- as.matrix(dist(X.match[in.sub,,drop = FALSE], method = k2k.method, p = mpower))

      #Put smaller group on rows
      d.rows <- which(rownames(dist.mat) %in% names(treat[in.sub])[treat[in.sub] == s])
      dist.mat <- dist.mat[d.rows, -d.rows, drop = FALSE]

      #For each member of group on row, find closest remaining pair from cols
      while (all(dim(dist.mat) > 0)) {
        extra.sub <- extra.sub + 1

        closest <- which.min(dist.mat[1,])
        subclass[c(rownames(dist.mat)[1], colnames(dist.mat)[closest])] <- extra.sub

        #Drop already paired units from dist.mat
        dist.mat <- dist.mat[-1,-closest, drop = FALSE]
      }

      #If any unmatched units remain, give them NA subclass
      if (any(dim(dist.mat) > 0)) subclass[unlist(dimnames(dist.mat))] <- NA_integer_

    }
  }

  subclass <- factor(subclass, nmax = extra.sub)


  names(subclass) <- names(treat)

  return(subclass)
}