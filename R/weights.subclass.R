weights.subclass <- function(psclass, treat, estimand = "ATT") {

  NAsub <- is.na(psclass)

  i1 <- treat == 1 & !NAsub
  i0 <- treat == 0 & !NAsub

  weights <- setNames(rep(0, length(treat)), names(treat))

  unique.sub <- unique(psclass[!NAsub], nmax = length(treat)/2)

  treated_by_sub <- vapply(unique.sub, function(s) sum(psclass[i1] == s), numeric(1L))
  control_by_sub <- vapply(unique.sub, function(s) sum(psclass[i0] == s), numeric(1L))
  total_by_sub <- treated_by_sub + control_by_sub

  psclass <- as.character(psclass)

  if (estimand == "ATT") {
    weights[i1] <- 1
    weights[i0] <- (treated_by_sub/control_by_sub)[psclass[i0]]

    #Weights average 1
    weights[i0] <- weights[i0]*sum(i0)/sum(weights[i0])
  }
  else if (estimand == "ATC") {
    weights[i1] <- (control_by_sub/treated_by_sub)[psclass[i1]]
    weights[i0] <- 1

    #Weights average 1
    weights[i1] <- weights[i1]*sum(i1)/sum(weights[i1])
  }
  else if (estimand == "ATE") {
    weights[i1] <- (total_by_sub/treated_by_sub)[psclass[i1]]
    weights[i0] <- (total_by_sub/control_by_sub)[psclass[i0]]

    #Weights average 1
    weights[i1] <- weights[i1]*sum(i1)/sum(weights[i1])
    weights[i0] <- weights[i0]*sum(i0)/sum(weights[i0])
  }

  if (sum(weights)==0)
    stop("No units were matched")
  else if (sum(weights[treat == 1])==0)
    stop("No treated units were matched")
  else if (sum(weights[treat == 0])==0)
    stop("No control units were matched")

  return(weights)
}
