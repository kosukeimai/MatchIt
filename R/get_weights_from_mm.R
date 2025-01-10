get_weights_from_mm <- function(match.matrix, treat, focal = NULL) {

  if (!is.integer(match.matrix)) {
    match.matrix <- charmm2nummm(match.matrix, treat)
  }

  weights <- weights_matrixC(match.matrix, treat, focal)

  if (all_equal_to(weights, 0)) {
    .err("no units were matched")
  }

  if (all_equal_to(weights[treat == 1], 0)) {
    .err("no treated units were matched")
  }

  if (all_equal_to(weights[treat == 0], 0)) {
    .err("no control units were matched")
  }

  setNames(weights, names(treat))
}