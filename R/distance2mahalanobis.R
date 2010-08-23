distance2mahalanobis <- function(formula, data, ...) {
  X <- model.matrix(formula, data)
  ## Placeholder where real work is done on a unit by unit basis
  distance <- rep(1, nrow(X))
  return(list(model = NULL,
              distance = distance))
}

