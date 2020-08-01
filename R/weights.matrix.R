weights.matrix <- function(match.matrix, treat) {

  lab1 <- rownames(match.matrix)
  lab0 <- names(treat)[treat == 0]

  weights <- setNames(rep(0, length(treat)), names(treat))

  NAmm <- is.na(match.matrix)
  num.matches <- rowSums(!NAmm)

  weights[lab1][num.matches[lab1] > 0] <- 1

  for (i in lab0[lab0 %in% match.matrix]) {
    in.match <- rowSums(!NAmm & match.matrix == i) > 0
    if (any(in.match)) weights[i] <- sum(1/num.matches[in.match])
  }

  if (sum(weights[lab0]) > 0) {
    weights[lab0] <- weights[lab0]*sum(weights[lab0] > 0)/sum(weights[lab0])
  }

  if (sum(weights)==0)
    stop("No units were matched")
  else if (sum(weights[lab1])==0)
    stop("No treated units were matched")
  else if (sum(weights[lab0])==0)
    stop("No control units were matched")

  return(weights)
}
