matchit.qqplot <- function(object, which.subclass = NULL,
                           interactive = TRUE, which.xs = NULL,...) {

  #Create covariate matrix; include exact and mahvars
  if (!is.null(which.xs)) {
    if (!is.vector(which.xs, "character")) {
      stop("which.xs should be a vector of variables for which balance is to be displayed.", call. = FALSE)
    }
    which.xs.f <- terms(reformulate(which.xs))

    X <- model.matrix(update(which.xs.f, NULL ~ . + 1), data = object$X,
                      contrasts.arg = lapply(Filter(is.factor, object$X[, all.vars(which.xs.f), drop = FALSE]),
                                             function(v) contrasts(v, contrasts = FALSE)))[,-1,drop = FALSE]
  }
  else {
    #Create covariate matrix; include exact and mahvars
    X <- model.matrix(update(object$formula, NULL ~ . + 1), data = object$X,
                      contrasts.arg = lapply(Filter(is.factor, object$X[, all.vars(object$formula), drop = FALSE]),
                                             function(v) contrasts(v, contrasts = FALSE)))[,-1,drop = FALSE]
    if (!is.null(object$exact)) {
      Xexact <- model.matrix(update(object$exact, NULL ~ . + 1), data = object$X,
                             contrasts.arg = lapply(Filter(is.factor, object$X[, all.vars(object$exact), drop = FALSE]),
                                                    function(v) contrasts(v, contrasts = FALSE)))[,-1,drop = FALSE]
      X <- cbind(X, Xexact[,setdiff(colnames(Xexact), colnames(X)), drop = FALSE])
    }

    if (!is.null(object$mahvars)) {
      Xmahvars <- model.matrix(update(object$mahvars, NULL ~ . + 1), data = object$X,
                               contrasts.arg = lapply(Filter(is.factor, object$X[, all.vars(object$mahvars), drop = FALSE]),
                                                      function(v) contrasts(v, contrasts = FALSE)))[,-1,drop = FALSE]
      X <- cbind(X, Xmahvars[,setdiff(colnames(Xmahvars), colnames(X)), drop = FALSE])
    }
  }

  discrete.cutoff <- 5

  t <- object$treat

  if (!is.null(object$subclass) && !is.null(which.subclass)) {
    if (!is.atomic(which.subclass) || !is.numeric(which.subclass) || length(which.subclass) > 1) {
      stop("The argument to subclass must be NULL or the index of the single subclass for which to display covariate distributions.", call. = FALSE)
    }
    if (!which.subclass %in% object$subclass[!is.na(object$subclass)]) {
      stop("The argument supplied to subclass is not the index of any subclass in the matchit object.", call. = FALSE)
    }
    w <- as.numeric(object$subclass == which.subclass & !is.na(object$subclass))
  }
  else {
    w <- object$weights
    if (is.null(w)) w <- rep(1, length(t))
  }

  for (i in 0:1) w[t == i] <- w[t==i]/sum(w[t==i])

  n1 <- sum(t == 1)
  n0 <- length(t) - n1

  wn1 <- sum(w[t==1] > 0)
  wn0 <- sum(w[t==0] > 0)

  varnames <- colnames(X)

  .pardefault <- par(no.readonly = TRUE)
  on.exit(par(.pardefault))

  oma <- c(2.25, 0, 3.75, 1.5)
  opar <- par(mfrow = c(3, 3), mar = rep.int(1/2, 4), oma = oma)

  for (i in seq_along(varnames)){
    x <- X[,i]

    plot(x, type= "n" , axes=FALSE)
    if(((i-1)%%3)==0){
      htext <- "QQ Plots"
      if (!is.null(which.subclass)){
        htext <- paste0(htext, paste0(" (Subclass ", which.subclass,")"))
      }
      mtext(htext, 3, 2, TRUE, 0.5, cex=1.1, font = 2)
      mtext("All", 3, .25, TRUE, 0.5, cex=1, font = 1)
      mtext("Matched", 3, .25, TRUE, 0.83, cex=1, font = 1)
      mtext("Control Units", 1, 0, TRUE, 2/3, cex=1, font = 1)
      mtext("Treated Units", 4, 0, TRUE, 0.5, cex=1, font = 1)
    }
    par(usr = c(0, 1, 0, 1))
    l.wid <- strwidth(varnames, "user")
    cex.labels <- max(0.75, min(1.45, 0.85/max(l.wid)))
    text(0.5, 0.5, varnames[i], cex = cex.labels)

    ord <- order(x)
    x_ord <- x[ord]
    t_ord <- t[ord]

    u <- unique(x_ord)

    #Need to interpolate larger group to be same size as smaller group

    #Unmatched sample
    x1 <- x_ord[t_ord == 1]
    x0 <- x_ord[t_ord != 1]

    if (n1 < n0) {
      if (length(u) <= discrete.cutoff) {
        x0probs <- vapply(u, function(u_) mean(x0 == u_), numeric(1L))
        x0cumprobs <- c(0, cumsum(x0probs)[-length(u)], 1)
        x0 <- u[findInterval(seq_len(n1)/n1, x0cumprobs, rightmost.closed = TRUE)]
      }
      else {
        x0 <- approx(seq_len(n0)/n0, y = x0, xout = seq_len(n1)/n1, rule = 2,
                     method = "constant", ties = "ordered")$y
      }
    }
    else {
      if (length(u) <= discrete.cutoff) {
        x1probs <- vapply(u, function(u_) mean(x1 == u_), numeric(1L))
        x1cumprobs <- c(0, cumsum(x1probs)[-length(u)], 1)
        x1 <- u[findInterval(seq_len(n0)/n0, x1cumprobs, rightmost.closed = TRUE)]
      }
      else {
        x1 <- approx(seq_len(n1)/n1, y = x1, xout = seq_len(n0)/n0, rule = 2,
                     method = "constant", ties = "ordered")$y
      }
    }

    if (length(u) <= discrete.cutoff) {
      md <- min(diff(u))
      x0 <- jitter(x0, amount = .07*md)
      x1 <- jitter(x1, amount = .07*md)
    }

    rr <- range(c(x0, x1))
    plot(x0, x1, xlab = "", ylab = "", xlim = rr, ylim = rr, axes = FALSE, ...)
    abline(a = 0, b = 1)
    abline(a = (rr[2]-rr[1])*0.1, b = 1, lty = 2)
    abline(a = -(rr[2]-rr[1])*0.1, b = 1, lty = 2)
    axis(2)
    box()

    #Matched sample
    w_ord <- w[ord]
    w1 <- w_ord[t_ord == 1]
    w0 <- w_ord[t_ord != 1]

    x1 <- x_ord[t_ord == 1][w1 > 0]
    x0 <- x_ord[t_ord != 1][w0 > 0]

    if (wn1 < wn0) {
      if (length(u) <= discrete.cutoff) {
        x0probs <- vapply(u, function(u_) weighted.mean(x0 == u_, w0[w0 > 0]), numeric(1L))
        x0cumprobs <- c(0, cumsum(x0probs)[-length(u)], 1)
        x0 <- u[findInterval(cumsum(w1[w1 > 0]), x0cumprobs, rightmost.closed = TRUE)]
      }
      else {
        x0 <- approx(cumsum(w0[w0 > 0]), y = x0, xout = cumsum(w1[w1 > 0]), rule = 2,
                     method = "constant", ties = "ordered")$y
      }
    }
    else {
      if (length(u) <= discrete.cutoff) {
        x1probs <- vapply(u, function(u_) weighted.mean(x1 == u_, w1[w1 > 0]), numeric(1L))
        x1cumprobs <- c(0, cumsum(x1probs)[-length(u)], 1)
        x1 <- u[findInterval(cumsum(w0[w0 > 0]), x1cumprobs, rightmost.closed = TRUE)]
      }
      else {
        x1 <- approx(cumsum(w1[w1 > 0]), y = x1, xout = cumsum(w0[w0 > 0]), rule = 2,
                     method = "constant", ties = "ordered")$y
      }
    }

    if (length(u) <= discrete.cutoff) {
      md <- min(diff(u))
      x0 <- jitter(x0, amount = .1*md)
      x1 <- jitter(x1, amount = .1*md)
    }

    plot(x0, x1, xlab = "", ylab = "", xlim = rr, ylim = rr, axes = FALSE, ...)
    abline(a = 0, b = 1)
    abline(a = (rr[2]-rr[1])*0.1, b = 1, lty = 2)
    abline(a = -(rr[2]-rr[1])*0.1, b = 1, lty = 2)
    box()

    devAskNewPage(ask = interactive)
  }
  devAskNewPage(ask = FALSE)
}
