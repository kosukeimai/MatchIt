get_weights_from_subclass <- function(psclass, treat, estimand = "ATT") {

  NAsub <- is.na(psclass)

  i1 <- which(treat == 1 & !NAsub)
  i0 <- which(treat == 0 & !NAsub)

  if (is_null(i1)) {
    if (is_null(i0)) {
      .err("No units were matched")
    }

    .err("No treated units were matched")
  }
  else if (is_null(i0)) {
    .err("No control units were matched")
  }

  weights <- setNames(rep(0.0, length(treat)), names(treat))

  if (!is.factor(psclass)) {
    psclass <- factor(psclass, nmax = min(length(i1), length(i0)))
  }

  treated_by_sub <- tabulate(psclass[i1], nlevels(psclass))
  control_by_sub <- tabulate(psclass[i0], nlevels(psclass))

  psclass <- unclass(psclass)

  if (estimand == "ATT") {
    weights[i1] <- 1.0
    weights[i0] <- (treated_by_sub/control_by_sub)[psclass[i0]]
  }
  else if (estimand == "ATC") {
    weights[i1] <- (control_by_sub/treated_by_sub)[psclass[i1]]
    weights[i0] <- 1.0
  }
  else if (estimand == "ATE") {
    weights[i1] <- 1.0 + (control_by_sub/treated_by_sub)[psclass[i1]]
    weights[i0] <- 1.0 + (treated_by_sub/control_by_sub)[psclass[i0]]
  }

  weights
}

# get_weights_from_subclass2 <- function(psclass, treat, estimand = "ATT") {
#
#   weights <- weights_subclassC(psclass, treat,
#                                switch(estimand, "ATT" = 1, "ATC" = 0, NULL))
#
#   if (sum(weights) == 0)
#     .err("No units were matched")
#   if (sum(weights[treat == 1]) == 0)
#     .err("No treated units were matched")
#   if (sum(weights[treat == 0]) == 0)
#     .err("No control units were matched")
#
#   weights
# }