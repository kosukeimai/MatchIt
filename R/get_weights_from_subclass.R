get_weights_from_subclass <- function(subclass, treat, estimand = "ATT") {

  NAsub <- is.na(subclass)

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

  w <- rep_with(0, treat)

  if (!is.factor(subclass)) {
    subclass <- factor(subclass, nmax = min(length(i1), length(i0)))
  }

  treated_by_sub <- tabulate(subclass[i1], nlevels(subclass))
  control_by_sub <- tabulate(subclass[i0], nlevels(subclass))

  subclass <- unclass(subclass)

  if (estimand == "ATT") {
    w[i1] <- 1
    w[i0] <- (treated_by_sub / control_by_sub)[subclass[i0]]
  }
  else if (estimand == "ATC") {
    w[i1] <- (control_by_sub / treated_by_sub)[subclass[i1]]
    w[i0] <- 1
  }
  else if (estimand == "ATE") {
    w[i1] <- 1 + (control_by_sub / treated_by_sub)[subclass[i1]]
    w[i0] <- 1 + (treated_by_sub / control_by_sub)[subclass[i0]]
  }

  w
}

# get_weights_from_subclass2 <- function(subclass, treat, estimand = "ATT") {
#
#   weights <- weights_subclassC(subclass, treat,
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