weights.matrix <- function(match.matrix, treat) {

  weights <- weights_matrixC(match.matrix, treat)

  if (sum(weights)==0)
    stop("No units were matched.", call. = FALSE)
  else if (sum(weights[treat == 1])==0)
    stop("No treated units were matched.", call. = FALSE)
  else if (sum(weights[treat == 0])==0)
    stop("No control units were matched.", call. = FALSE)

  return(setNames(weights, names(treat)))
}