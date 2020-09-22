nn_match <- function(treat, ord, ratio = 1, replace = FALSE, discarded, distance = NULL, ex = NULL, caliper.dist = NULL,
                    caliper.covs = NULL, caliper.covs.mat = NULL, mahcovs = NULL, mahSigma_inv = NULL) {

  n1 <- sum(treat == 1)

  ind <- seq_along(treat)
  ind1 <- ind[treat == 1]
  ind0 <- ind[treat == 0]

  if (!is.null(ex)) {
    ex1 <- ex[treat == 1]
    ex0 <- ex[treat == 0]
  }

  if (!is.null(distance)) {
    d1 <- distance[treat == 1]
    d0 <- distance[treat == 0]
  }

  max.ratio <- max(ratio)

  mm <- matrix(NA_integer_, nrow = n1, ncol = max.ratio)

  matched <- discarded

  ord_ <- ord[!discarded[ord]] #Only non-discarded

  for (r in seq_len(max.ratio)) {
    for (ord_i in ord_[ratio[ord_] >= r]) {

      c.eligible <- !matched & treat == 0

      if (!replace) {
        if (!any(c.eligible)) break
      }
      else if (r > 1) {
        #If replace = T and r > 1, don't rematch to same control unit
        c.eligible[mm[ord_i, seq_len(r - 1)]] <- FALSE
      }

      if (!any(c.eligible)) next

      if (!is.null(ex)) {
        c.eligible[c.eligible][ex0[c.eligible] != ex1[ord_i]] <- FALSE
      }

      if (!any(c.eligible)) next

      #Get distances among eligible and apply caliper

      ps.diff <- NULL

      #PS caliper
      if (length(caliper.dist) > 0) {
        ps.diff <- abs(d1[ord_i] - distance[c.eligible])
        c.eligible[c.eligible][ps.diff > caliper.dist] <- FALSE
      }

      if (!any(c.eligible)) next

      #Covariate caliper
      if (length(caliper.covs) > 0) {
        for (x in names(caliper.covs)) {
          calcov.diff <- abs(caliper.covs.mat[ind1[ord_i], x] - caliper.covs.mat[c.eligible, x])
          c.eligible[c.eligible][calcov.diff > caliper.covs[x]] <- FALSE
          if (!any(c.eligible)) break
        }
      }

      if (!any(c.eligible)) next

      if (length(mahcovs) == 0) {
        #PS matching
        distances <- if (is.null(ps.diff)) abs(d1[ord_i] - distance[c.eligible]) else ps.diff
      }
      else {
        #MD matching
        distances <- sqrt(mahalanobis(mahcovs[c.eligible, ,drop = FALSE],
                                      mahcovs[ind1[ord_i],],
                                      cov = mahSigma_inv, inverted = TRUE))
      }

      #Assign match
      ##Resolve ties by selecting the first unit
      mm[ord_i, r] <- which(c.eligible)[which.min(distances)]

      if (!replace) matched[mm[ord_i, r]] <- TRUE

    }
  }

  return(mm)
}