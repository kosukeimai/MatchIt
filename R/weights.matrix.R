weights.matrix <- function(match.matrix, treat) {

  if (!is.integer(match.matrix)) match.matrix <- charmm2nummm(match.matrix, treat)

  weights <- weights_matrixC(match.matrix, treat)

  if (sum(weights) == 0)
    .err("No units were matched")
  if (sum(weights[treat == 1]) == 0)
    .err("No treated units were matched")
  if (sum(weights[treat == 0]) == 0)
    .err("No control units were matched")

  setNames(weights, names(treat))
}