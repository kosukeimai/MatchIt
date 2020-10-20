weights.matrix <- function(match.matrix, treat) {

  in_t1 <- which(treat == 1)
  in_t0 <- which(treat != 1)

  weights <- rep(0, length(treat))

  match.matrix[is.na(match.matrix)] <- 0L
  num.matches <- rowSums(match.matrix != 0)

  weights[in_t1[num.matches > 0]] <- 1

  in_t0_matched <- in_t0[in_t0 %in% match.matrix]

  for (i in in_t0_matched) {
    in.match <- which(rowSums(match.matrix == i) > 0)
    weights[i] <- sum(1/num.matches[in.match])
  }

  if (sum(weights[in_t0]) > 0) {
    weights[in_t0] <- weights[in_t0]*sum(weights[in_t0] > 0)/sum(weights[in_t0])
  }

  if (sum(weights)==0)
    stop("No units were matched.", call. = FALSE)
  else if (sum(weights[in_t1])==0)
    stop("No treated units were matched.", call. = FALSE)
  else if (sum(weights[in_t0])==0)
    stop("No control units were matched.", call. = FALSE)

  return(setNames(weights, names(treat)))
}