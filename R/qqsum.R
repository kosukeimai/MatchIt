## Function for QQ summary stats
qqsum <- function(x, t, w = NULL, standardize = FALSE) {
  #x = variable, t = treat, w = weights

  n.obs <- length(x)

  if (is.null(w)) w <- rep(1, n.obs)

  if (all(x == 0 | x == 1)) {
    #For binary variables, just difference in means
    ediff <- abs(weighted.mean(x[t == t[1]], w[t == t[1]]) - weighted.mean(x[t != t[1]], w[t != t[1]]))
    return(c(meandiff = ediff, meddiff = ediff, maxdiff = ediff))
  }
  else {
    for (i in unique(t, nmax = 2)) w[t==i] <- w[t==i]/sum(w[t==i])

    ord <- order(x)
    x_ord <- x[ord]
    w_ord <- w[ord]
    t_ord <- t[ord]

    if (standardize) {
      #Difference between ecdf of x for each group
      w_ord_ <- w_ord
      w_ord_[t_ord==t_ord[1]] <- -w_ord_[t_ord==t_ord[1]]
      ediff <- abs(cumsum(w_ord_))[c(diff(x_ord) != 0, TRUE)]
    }
    else {
      #Horizontal distance of ecdf between groups
      #Need to interpolate larger group to be same size as smaller group

      u <- unique(x_ord)

      wn1 <- sum(w[t == t_ord[1]] > 0)
      wn0 <- sum(w[t != t_ord[1]] > 0)

      w1 <- w_ord[t_ord == t_ord[1]]
      w0 <- w_ord[t_ord != t_ord[1]]

      x1 <- x_ord[t_ord == t_ord[1]][w1 > 0]
      x0 <- x_ord[t_ord != t_ord[1]][w0 > 0]

      if (wn1 < wn0) {
        if (length(u) <= 5) {
          x0probs <- vapply(u, function(u_) weighted.mean(x0 == u_, w0[w0 > 0]), numeric(1L))
          x0cumprobs <- c(0, cumsum(x0probs)[-length(u)], 1)
          x0 <- u[findInterval(cumsum(w1[w1 > 0]), x0cumprobs, rightmost.closed = TRUE)]
        }
        else {
          x0 <- approx(cumsum(w0[w0 > 0]), y = x0, xout = cumsum(w1[w1 > 0]), rule = 2,
                       method = "constant", ties = "ordered")$y
        }
      }
      else if (wn1 > wn0) {
        if (length(u) <= 5) {
          x1probs <- vapply(u, function(u_) weighted.mean(x1 == u_, w1[w1 > 0]), numeric(1L))
          x1cumprobs <- c(0, cumsum(x1probs)[-length(u)], 1)
          x1 <- u[findInterval(cumsum(w0[w0 > 0]), x1cumprobs, rightmost.closed = TRUE)]
        }
        else {
          x1 <- approx(cumsum(w1[w1 > 0]), y = x1, xout = cumsum(w0[w0 > 0]), rule = 2,
                       method = "constant", ties = "ordered")$y
        }
      }
      ediff <- abs(x1 - x0)
    }
    return(c(meandiff = mean(ediff), meddiff = median(ediff), maxdiff = max(ediff)))
  }
}