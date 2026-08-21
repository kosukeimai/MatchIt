get_weights_from_subclass <- function(subclass, treat, estimand = "ATT", s.weights = NULL) {

  if (is_null(s.weights)) {
    s.weights <- rep_with(1, treat)
  }

  NAsub <- is.na(subclass)

  i1 <- treat == 1 & !NAsub
  i0 <- treat == 0 & !NAsub

  if (sum(s.weights[i1]) == 0) {
    if (sum(s.weights[i0]) == 0) {
      arg::err("no units were matched")
    }

    arg::err("no treated units were matched")
  }
  else if (sum(s.weights[i0]) == 0) {
    arg::err("no control units were matched")
  }

  w <- rep_with(0, treat)

  if (!is.factor(subclass)) {
    subclass <- factor(subclass, nmax = min(sum(i1), sum(i0)))
  }

  subclass_mass <- matrix(0, nrow = 2L, ncol = nlevels(subclass),
                          dimnames = list(NULL, levels(subclass)))

  for (s in levels(subclass)) {
    in_s <- subclass == s
    subclass_mass[, s] <- c(sum(s.weights[in_s & i0]),
                            sum(s.weights[in_s & i1]))
  }

  subclass <- levels(subclass)[unclass(subclass)]

  if (estimand == "ATT") {
    w[i1] <- 1
    w[i0] <- (subclass_mass[2L, ] / subclass_mass[1L, ])[subclass[i0]]
  }
  else if (estimand == "ATC") {
    w[i1] <- (subclass_mass[1L, ] / subclass_mass[2L, ])[subclass[i1]]
    w[i0] <- 1
  }
  else if (estimand == "ATE") {
    w[i1] <- 1 + (subclass_mass[1L, ] / subclass_mass[2L, ])[subclass[i1]]
    w[i0] <- 1 + (subclass_mass[2L, ] / subclass_mass[1L, ])[subclass[i0]]
  }

  w
}