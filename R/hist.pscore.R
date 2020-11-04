hist.pscore <- function(x, xlab="Propensity Score", freq = FALSE, ...){
  .pardefault <- par(no.readonly = TRUE)
  on.exit(par(.pardefault))

  treat <- x$treat
  pscore <- x$distance[!is.na(x$distance)]
  s.weights <- if (is.null(x$s.weights)) rep(1, length(treat)) else x$s.weights
  weights <- x$weights * s.weights
  matched <- weights != 0
  q.cut <- x$q.cut

  minp <- min(pscore)
  maxp <- max(pscore)
  ratio <- x$call$ratio

  if (is.null(ratio)) ratio <- 1

  for (i in unique(treat, nmax = 2)) {
    if (freq) s.weights[treat == i] <- s.weights[treat == i]/mean(s.weights[treat == i])
    else s.weights[treat == i] <- s.weights[treat == i]/sum(s.weights[treat == i])

    if (freq) weights[treat == i] <- weights[treat == i]/mean(weights[treat == i])
    else weights[treat == i] <- weights[treat == i]/sum(weights[treat == i])
  }

  ylab <- if (freq) "Count" else "Proportion"

  par(mfrow = c(2,2))
  # breaks <- pretty(na.omit(pscore), 10)
  breaks <- seq(minp, maxp, length = 11)
  xlim <- range(breaks)

  for (n in c("Raw Treated", "Matched Treated", "Raw Control", "Matched Control")) {
    if (startsWith(n, "Raw")) w <- s.weights
    else w <- weights

    if (endsWith(n, "Treated")) t <- 1
    else t <- 0

    #Create histogram using weights
    #Manually assign density, which is used as height of the bars. The scaling
    #of the weights above determine whether they are "counts" or "proportions".
    #Regardless, set freq = FALSE in plot() to ensure density is used for bar
    #height rather than count.
    pm <- hist(pscore[treat==t], plot = FALSE, breaks = breaks)
    pm[["density"]] <- vapply(seq_len(length(pm$breaks) - 1), function(i) {
      sum(w[treat == t & pscore >= pm$breaks[i] & pscore < pm$breaks[i+1]])
    }, numeric(1L))
    plot(pm, xlim = xlim, xlab = xlab, main = n, ylab = ylab,
         freq = FALSE, col = "lightgray", ...)

    if (!startsWith(n, "Raw") && !is.null(q.cut)) abline(v = q.cut, lty=2)
  }

}
