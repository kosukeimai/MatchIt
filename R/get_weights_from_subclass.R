get_weights_from_subclass <- function(psclass, treat, estimand = "ATT") {

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
    weights[i0] <- .make_sum_to_n(weights[i0])
  }
  else if (estimand == "ATC") {
    weights[i1] <- (control_by_sub/treated_by_sub)[psclass[i1]]
    weights[i0] <- 1

    #Weights average 1
    weights[i1] <- .make_sum_to_n(weights[i1])
  }
  else if (estimand == "ATE") {
    weights[i1] <- (total_by_sub/treated_by_sub)[psclass[i1]]
    weights[i0] <- (total_by_sub/control_by_sub)[psclass[i0]]

    #Weights average 1
    weights[i1] <- .make_sum_to_n(weights[i1])
    weights[i0] <- .make_sum_to_n(weights[i0])
  }

  if (sum(weights) == 0)
    .err("No units were matched")
  if (sum(weights[treat == 1]) == 0)
    .err("No treated units were matched")
  if (sum(weights[treat == 0]) == 0)
    .err("No control units were matched")

  weights
}
