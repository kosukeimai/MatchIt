get_weights_from_mm <- function(match.matrix, treat, focal = NULL) {

  if (!is.integer(match.matrix)) {
    match.matrix <- charmm2nummm(match.matrix, treat)
  }

  w <- weights_matrixC(match.matrix, treat, focal)

  if (all_equal_to(w, 0)) {
    .err("no units were matched")
  }

  if (all_equal_to(w[treat == 1], 0)) {
    .err("no treated units were matched")
  }

  if (all_equal_to(w[treat == 0], 0)) {
    .err("no control units were matched")
  }

  setNames(w, names(treat))
}