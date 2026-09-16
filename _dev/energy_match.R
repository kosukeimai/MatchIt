energy_match <- function(treat, dist_mat, estimand = "ATT", ratio = NULL, s.weights = NULL, verbose = FALSE) {

  if (estimand == "ATT") {
    out <- energy_match_att(treat, dist_mat, verbose)
  }
  else if (estimand == "ATC") {
    out <- energy_match_att(1 - treat, dist_mat, verbose)
  }
  else if (estimand == "ATE") {
    out <- energy_match_ate(treat, dist_mat, NULL, verbose)
  }

  drop_order <- out[[1]]
  edists <- out[[2]]

  last_drop <- which.min(edists)

  weights <- rep(1, length(treat))
  weights[drop_order[seq_len(last_drop)[-1]]] <- 0

  return(weights)
}

pow <- function(x, n) x^n

energy_match_att <- function(t, d, verbose = FALSE) {
  treated_ind = which(t == 1)
  control_ind = which(t == 0)

  n1 <- length(treated_ind)
  n0 <- n0_ <- length(control_ind)

  drop_order <- rep(NA_integer_, n0 + 1)
  Ys <- rep(NA_real_, n0 + 1)

  d10 <- d[treated_ind, control_ind, drop = FALSE]
  cSd10 <- colSums(d10)
  Sd10 <- sum(cSd10)

  d11 <- d[treated_ind, treated_ind, drop = FALSE]
  Sd11 <- sum(d11)

  d00 <- d[control_ind, control_ind, drop = FALSE]
  cSd00 <- colSums(d00)
  Sd00 <- sum(cSd00)
  d00i <- seq_len(n0)

  Ys[1] <- 2*Sd10/(n1*n0) - Sd11/(n1*n1) - Sd00/(n0*n0)

  min_edist <- Ys[1]
  edist <- numeric(1)

  control_contributions <-
    (2/(n0*n1) - 2/((n0-1)*n1))*Sd10 +
    (-1/(n0^2) + 1/((n0-1)^2))*Sd00 +
    (2/((n0-1)*n1))*cSd10 -
    (2/((n0-1)^2))*cSd00

  k <- 2
  repeat {
    # Find which units have the largest contribution to the energy distance
    drop_control <- which.max(control_contributions)
    control_ind_to_drop <- control_ind[drop_control]

    control_ind <- control_ind[-drop_control]

    n0 <- length(control_ind)

    # //If removed units are last units, don't remove and stop
    if (n0 - 1 <= 0) break

    # // p.update(n0_ - n0)

    # //Compute new edist with dropped units removed

    d00i_control_dropped <- d00i[drop_control]
    d00i <- d00i[-drop_control]

    # //d10 contributions
    Sd10 <- Sd10 - sum(d10[,d00i_control_dropped, drop = FALSE])
    cSd10 <- cSd10[-drop_control]

    # //d00 contributions
    Sd00 <- Sd00 - 2*sum(d00[d00i, d00i_control_dropped, drop = FALSE]) -
      sum(d00[d00i_control_dropped, d00i_control_dropped, drop = FALSE])

    cSd00 <- cSd00[-drop_control]
    cSd00 <- cSd00 - colSums(d00[d00i_control_dropped, d00i, drop = FALSE])

    edist <- 2*Sd10/(n1*n0) - Sd11/(n1*n1) - Sd00/(n0*n0)

    # //After passing sample size threshold of 90% of original N, stop if
    # //new edist is larger than smallest edist
    if (edist < min_edist) min_edist <- edist
    else if (n0/n0_ < .9 && edist - min_edist > .2*(Ys[1] - min_edist)) break

    # //Record new edist and units dropped
    Ys[k] <- edist
    drop_order[k] <- control_ind_to_drop

    # //Compute unit contributions after having dropped units
    control_contributions <-
      (2/(n0*n1) - 2/((n0-1)*n1))*Sd10 +
      (-1/(n0^2) + 1/((n0-1)^2))*Sd00 +
      (2/((n0-1)*n1))*cSd10 -
      (2/((n0-1)^2))*cSd00

    k <- k + 1
  }

  list(drop_order = drop_order, edist = Ys)
}

energy_match_ate <- function(t, d, ratio = NULL, verbose = FALSE) {
  #Energy distance formula:
  #2/n1n0 sum(D[t,c]) - 1/n1n1 sum(D[t,t]) - 1/n0n0 sum(D[c,c])

  n <- length(t)

  treated_ind = which(t == 1)
  control_ind = which(t == 0)

  n1 <- length(treated_ind)
  n0 <- length(control_ind)

  drop_order <- rep(NA_integer_, n + 1)
  Ys <- rep(NA_real_, n + 1)

  #Compute each unit's contribution to energy dist
  #Find unit that contributes maximally to energy distance
  #Remove that unit

  # if (verbose) {
  #   pb <- pbapply::startpb(min = 0, max = N)
  # }

  Sdff <- sum(d)

  df1 <- d[,treated_ind, drop = FALSE]
  cSdf1 <- colSums(df1)
  Sdf1 <- sum(cSdf1)

  df0 <- d[,control_ind, drop = FALSE]
  cSdf0 <- colSums(df0)
  Sdf0 <- Sdff - Sdf1

  d10 <- df0[treated_ind,, drop = FALSE]
  rSd10 <- rowSums(d10)
  cSd10 <- colSums(d10)
  Sd10 <- sum(rSd10)

  d11 <- df1[treated_ind,,drop = FALSE]
  cSd11 <- colSums(d11)
  Sd11 <- sum(cSd11)
  d11i <- seq_len(ncol(d11))

  d00 <- df0[control_ind,, drop = FALSE]
  cSd00 <- colSums(d00)
  Sd00 <- sum(cSd00)
  d00i <- seq_len(ncol(d00))

  Ys[1] <- 2*Sd10/(n1*n0) - Sd11/(n1*n1) - Sd00/(n0*n0) +
    2*Sdf1/(n1*n) - Sd11/(n1*n1) - Sdff/(n*n) +
    2*Sdf0/(n0*n) - Sd00/(n0*n0) - Sdff/(n*n)

  min_edist <- Ys[1]

  #Difference between edist and edist with each unit
  treated_contributions <- (2/(n1*n0) - 2/((n1-1)*n0))*Sd10 +
    (-1/(n1^2) + 1/((n1-1)^2))*Sd11 +
    (2/((n1-1)*n0))*rSd10 -
    (2/((n1-1)^2))*cSd11 +
    (2/(n1*n) - 2/((n1-1)*n))*Sdf1 +
    (-1/(n1^2) + 1/((n1-1)^2))*Sd11 +
    (2/((n1-1)*n))*cSdf1 -
    (2/((n1-1)^2))*cSd11
  control_contributions <- (2/(n0*n1) - 2/((n0-1)*n1))*Sd10 +
    (-1/(n0^2) + 1/((n0-1)^2))*Sd00 +
    (2/((n0-1)*n1))*cSd10 -
    (2/((n0-1)^2))*cSd00 +
    (2/(n0*n) - 2/((n0-1)*n))*Sdf0 +
    (-1/(n0^2) + 1/((n0-1)^2))*Sd00 +
    (2/((n0-1)*n))*cSdf0 -
    (2/((n0-1)^2))*cSd00

  k <- 2
  repeat {

    #Find which units have the largest contribution to the energy distance
    if (!is.null(ratio)) {
      if (n0 <= n1 * ratio) {
        largest_contribution <- max(treated_contributions)
        drop_treated <- which.min(abs(treated_contributions - largest_contribution))
        drop_control <- integer(0L)
      }
      else {
        largest_contribution <- max(control_contributions)
        drop_treated <- integer(0L)
        drop_control <- which.min(abs(control_contributions - largest_contribution))
      }
    }
    else {
      largest_contribution <- max(max(treated_contributions), max(control_contributions))

      drop_treated <- which(treated_contributions == largest_contribution)
      drop_control <- which(control_contributions == largest_contribution)

      if (length(drop_treated) > 0 && length(drop_control) > 0) {
        if (length(treated_contributions) > length(control_contributions)) {
          drop_control <- integer(0L)
          drop_treated <- drop_treated[1]
        }
        else {
          drop_control <- drop_control[1]
          drop_treated <- integer(0L)
        }
      }
      else if (length(drop_treated) > 0) {
        drop_treated <- drop_treated[1]
      }
      else {
        drop_control <- drop_control[1]
      }
    }

    treated_ind_to_drop <- treated_ind[drop_treated]
    control_ind_to_drop <- control_ind[drop_control]

    #Remove those units and decrease
    if (length(drop_treated) > 0) {
      treated_ind <- treated_ind[-drop_treated]
      n1 <- length(treated_ind)
    }
    if (length(drop_control) > 0) {
      control_ind <- control_ind[-drop_control]
      n0 <- length(control_ind)
    }

    #If removed units are last units, don't remove and stop
    if (n1-1 <= 0 || n0-1 <= 0) break

    # if (verbose) {
    #   pbapply::setpb(pb, n - (n0 + n1))
    # }

    #Compute new edist with dropped units removed
    if (length(drop_treated) > 0) {
      d11i_treated_dropped <- d11i[drop_treated]
      d11i_treated_kept <- d11i[-drop_treated]
    }
    else {
      d11i_treated_dropped <- integer(0L)
      d11i_treated_kept <- d11i
    }
    if (length(drop_control) > 0) {
      d00i_control_dropped <- d00i[drop_control]
      d00i_control_kept <- d00i[-drop_control]
    }
    else {
      d00i_control_dropped <- integer(0L)
      d00i_control_kept <- d00i
    }

    #d10 contributions
    if (length(drop_treated) > 0 && length(drop_control) > 0) {
      Sd10 <- Sd10 - sum(d10[d11i_treated_dropped, d00i_control_dropped]) -
        sum(d10[d11i_treated_dropped, d00i_control_kept]) -
        sum(d10[d11i_treated_kept, d00i_control_dropped])

      rSd10 <- rSd10[-drop_treated] - rowSums(d10[d11i_treated_kept, d00i_control_dropped, drop = FALSE])
      cSd10 <- cSd10[-drop_control] - colSums(d10[d11i_treated_dropped, d00i_control_kept, drop = FALSE])
    }
    else if (length(drop_treated) > 0) {
      Sd10 <- Sd10 - sum(d10[d11i_treated_dropped, d00i_control_kept])

      rSd10 <- rSd10[-drop_treated]
      cSd10 <- cSd10 - colSums(d10[d11i_treated_dropped, d00i_control_kept, drop = FALSE])
    }
    else if (length(drop_control) > 0) {
      Sd10 <- Sd10 - sum(d10[d11i_treated_kept, d00i_control_dropped])
      rSd10 <- rSd10 - rowSums(d10[d11i_treated_kept, d00i_control_dropped, drop = FALSE])
      cSd10 <- cSd10[-drop_control]
    }

    #d11 and df1 contributions
    if (length(drop_treated) > 0) {
      Sd11 <- Sd11 - 2*sum(d11[d11i_treated_kept, d11i_treated_dropped]) -
        sum(d11[d11i_treated_dropped, d11i_treated_dropped])
      cSd11 <- cSd11[-drop_treated] - colSums(d11[d11i_treated_dropped, d11i_treated_kept, drop = FALSE])
      d11i <- d11i_treated_kept

      Sdf1 <- Sdf1 - sum(df1[, d11i_treated_dropped])
      cSdf1 <- cSdf1[-drop_treated]
    }

    #d00 and df0 contributions
    if (length(drop_control) > 0) {
      Sd00 <- Sd00 - 2*sum(d00[d00i_control_kept, d00i_control_dropped]) -
        sum(d00[d00i_control_dropped, d00i_control_dropped])
      cSd00 <- cSd00[-drop_control] - colSums(d00[d00i_control_dropped, d00i_control_kept, drop = FALSE])
      d00i <- d00i_control_kept

      Sdf0 <- Sdf0 - sum(df0[, d00i_control_dropped])
      cSdf0 <- cSdf0[-drop_control]
    }

    edist <- 2*Sd10/(n1*n0) - Sd11/(n1*n1) - Sd00/(n0*n0) +
      2*Sdf1/(n1*n) - Sd11/(n1*n1) - Sdff/(n*n) +
      2*Sdf0/(n0*n) - Sd00/(n0*n0) - Sdff/(n*n)

    #After passing sample size threshold of 90% of original n, stop if
    #new edist is larger than smallest edist
    if (edist < min_edist) min_edist <- edist
    else if ((n1+n0)/n < .9 && edist-min_edist > .2*(Ys[1]-min_edist)) break

    #Record new edist and units dropped
    Ys[k] <- edist
    drop_order[k] <- c(treated_ind_to_drop, control_ind_to_drop)

    treated_contributions <- (2/(n1*n0) - 2/((n1-1)*n0))*Sd10 +
      (-1/(n1^2) + 1/((n1-1)^2))*Sd11 +
      (2/((n1-1)*n0))*rSd10 -
      (2/((n1-1)^2))*cSd11 +
      (2/(n1*n) - 2/((n1-1)*n))*Sdf1 +
      (-1/(n1^2) + 1/((n1-1)^2))*Sd11 +
      (2/((n1-1)*n))*cSdf1 -
      (2/((n1-1)^2))*cSd11
    control_contributions <- (2/(n0*n1) - 2/((n0-1)*n1))*Sd10 +
      (-1/(n0^2) + 1/((n0-1)^2))*Sd00 +
      (2/((n0-1)*n1))*cSd10 -
      (2/((n0-1)^2))*cSd00 +
      (2/(n0*n) - 2/((n0-1)*n))*Sdf0 +
      (-1/(n0^2) + 1/((n0-1)^2))*Sd00 +
      (2/((n0-1)*n))*cSdf0 -
      (2/((n0-1)^2))*cSd00

    k <- k + 1
  }

  Ys <- Ys[seq_len(k-1)]
  drop_order <- drop_order[seq_len(k-1)]

#   if (verbose) {
#     pbapply::setpb(pb, n)
#     pbapply::closepb(pb)
#   }

  return(list(drop_order = drop_order, edist = Ys))
}
