weights.subclass <- function(psclass, treat, estimand = "ATT") {

  NAsub <- is.na(psclass)

  i1 <- treat == 1 & !NAsub
  i0 <- treat == 0 & !NAsub

  weights <- setNames(rep(0, length(treat)), names(treat))

  if (!is.factor(psclass)) {
    psclass <- factor(psclass, nmax = min(sum(i1), sum(i0)))
    levels(psclass) <- seq_len(nlevels(psclass))
  }

  treated_by_sub <- setNames(tabulateC(psclass[i1], nlevels(psclass)), levels(psclass))
  control_by_sub <- setNames(tabulateC(psclass[i0], nlevels(psclass)), levels(psclass))

  total_by_sub <- treated_by_sub + control_by_sub

  psclass <- as.character(psclass)

  if (estimand == "ATT") {
    weights[i1] <- 1
    weights[i0] <- (treated_by_sub/control_by_sub)[psclass[i0]]

    #Weights average 1
    weights[i0] <- weights[i0]*sum(i0)/sum(weights[i0])
  }
  else if (estimand == "ATC") {
    weights[i1] <- (control_by_sub/treated_by_sub)[psclass[i1]]
    weights[i0] <- 1

    #Weights average 1
    weights[i1] <- weights[i1]*sum(i1)/sum(weights[i1])
  }
  else if (estimand == "ATE") {
    weights[i1] <- (total_by_sub/treated_by_sub)[psclass[i1]]
    weights[i0] <- (total_by_sub/control_by_sub)[psclass[i0]]

    #Weights average 1
    weights[i1] <- weights[i1]*sum(i1)/sum(weights[i1])
    weights[i0] <- weights[i0]*sum(i0)/sum(weights[i0])
  }

  if (sum(weights)==0)
    stop("No units were matched.", call. = FALSE)
  else if (sum(weights[treat == 1])==0)
    stop("No treated units were matched.", call. = FALSE)
  else if (sum(weights[treat == 0])==0)
    stop("No control units were matched.", call. = FALSE)

  return(weights)
}

#Not used yet, will be for multi-category treats
weights.subclass.multi <- function(psclass, treat, focal = NULL) {

  NAsub <- is.na(psclass)

  unique_treat <- sort(unique(treat))
  ind <- vapply(unique_treat, function(t) {
    treat == t & !NAsub
  }, logical(length(treat)))

  weights <- setNames(rep(0, length(treat)), names(treat))

  if (!is.factor(psclass)) {
    psclass <- factor(psclass, nmax = min(colSums(ind)))
    levels(psclass) <- seq_len(nlevels(psclass))
  }

  count_by_sub <- vapply(seq_along(unique_treat), function(i) {
    tabulateC(psclass[ind[,i]], nlevels(psclass))
  }, numeric(nlevels(psclass)))
  rownames(count_by_sub) <- levels(psclass)

  total_by_sub <- rowSums(count_by_sub)

  psclass <- as.character(psclass)

  if (is.null(focal)) {
    for (i in seq_along(unique_treat)) {
      psclass_i <- psclass[ind[,i]]
      weights[ind[,i]] <- total_by_sub[psclass_i]/count_by_sub[psclass_i,i]

      #Weights average 1
      weights[ind[,i]] <- weights[ind[,i]]*sum(ind[,i])/sum(weights[ind[,i]])
    }
  }
  else {
    focal_treat <- match(focal, unique_treat)
    weights[ind[,focal_treat]] <- 1

    for (i in seq_along(unique_treat)[-focal_treat]) {
      psclass_i <- psclass[ind[,i]]
      weights[ind[,i]] <- count_by_sub[psclass_i,focal_treat]/count_by_sub[psclass_i,i]

      #Weights average 1
      weights[ind[,i]] <- weights[ind[,i]]*sum(ind[,i])/sum(weights[ind[,i]])
    }
  }

  if (sum(weights)==0)
    stop("No units were matched.", call. = FALSE)
  else if (sum(weights[treat == 1])==0)
    stop("No treated units were matched.", call. = FALSE)
  else if (sum(weights[treat == 0])==0)
    stop("No control units were matched.", call. = FALSE)

  return(weights)
}
