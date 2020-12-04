cem_matchit <- function(treat, X, cutpoints = NULL, grouping = NULL, k2k = FALSE, k2k.method = "mahalanobis", mpower = 2, estimand = "ATT") {
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

  #Process cutpoints
  for (i in names(X)[sapply(X, is.numeric)]) {
    if (is.null(cutpoints) || !i %in% names(cutpoints)) breaks <- "sturges"
    else breaks <- cutpoints[[i]]

    if (is.character(breaks)) {
      breaks <- match_arg(tolower(breaks), c("sturges", "fd", "scott", "ss"))
      breaks <- switch(breaks,
                       sturges = nclass.Sturges(X[[i]]),
                       fd = nclass.FD(X[[i]]),
                       scott = nclass.scott(X[[i]]),
                       stop("unknown 'breaks' algorithm"))
    }
    if (is.numeric(breaks)) {
      if (length(breaks) == 1) {
        #cutpoints is number of bins, unlike in cem
        breaks <- seq(min(X[[i]]), max(X[[i]]), length = breaks + 1)
      }
      else {
        breaks <- sort(unique(breaks))
      }
    }
    X[[i]] <- findInterval(X[[i]], breaks, all.inside = TRUE)
  }

  #Exact match
  xx <- exactify(X, names(treat))
  cc <- intersect(xx[treat==1], xx[treat==0])

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