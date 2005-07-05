matchit2optimal <- function(treat, X, data, pscore, discarded,
                            ratio = 1, ...) {
  require(optmatch)

  ## optimal matching for undiscarded units
  ttt <- treat[!discarded]
  n0 <- length(ttt[ttt==0])
  n1 <- length(ttt[ttt==1])
  d1 <- pscore[ttt==1]
  d0 <- pscore[ttt==0]
  d <- matrix(0, ncol=n0, nrow=n1)
  tlabels <- rownames(d) <- names(ttt[ttt==1])
  clabels <- colnames(d) <- names(ttt[ttt==0])
  for (i in 1:n1) 
    d[i,] <- abs(d1[i]-d0)
  full <- fullmatch(d, min.controls = ratio,
                    max.controls = ratio, omit.fraction =
                    (n0-ratio*n1)/n0, ...)
  psclass <- full[pmatch(names(ttt), names(full))]
  psclass <- as.numeric(as.factor(psclass))
  names(psclass) <- names(ttt)

  mm <- matrix(0, nrow = n1, ncol = ratio, dimnames = list(tlabels, 1:ratio))
  for (i in 1:n1)
    mm[i,] <- names(which(psclass[tlabels[i]] == psclass[-pmatch(tlabels[i],
                                   names(psclass))]))

  if (any(discarded)) {
    ## add psclass = NA for discarded units
    tmp <- rep(NA, sum(discarded))
    names(tmp) <- names(treat[discarded])
    psclass <- c(psclass, tmp)[names(treat)]

    ## add match.matrix = NA for discarded units
    tdisc <- discarded[treat==1]
    if (any(tdisc)) {
      tmp <- matrix(NA, nrow = sum(tdisc), ncol= ratio,
                    dimnames = list(names(treat[treat==1 & discarded]),
                      1:ratio))
      mm <- rbind(mm, tmp)[names(treat[treat==1]),]
    }
  }

  ## calculate weights and return the results
  psweights <- weights.matrix(match.matrix = mm, treat=treat)$psweights
  res <- list(match.matrix = mm, psclass=psclass, psweights=psweights)

  class(res) <- "matchit"
  return(res)
}
