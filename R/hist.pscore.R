hist.pscore <- function(x, xlab="Propensity Score", freq = FALSE, ...){
  .pardefault <- par(no.readonly = TRUE)
  on.exit(par(.pardefault))

  treat <- x$treat
  pscore <- x$distance[!is.na(x$distance)]
  weights <- x$weights
  matched <- weights!=0
  q.cut <- x$q.cut
  cwt <- sqrt(weights)
  minp <- min(pscore)
  maxp <- max(pscore)
  ratio <- x$call$ratio

  if (is.null(ratio)) {ratio <- 1}

  for (i in unique(treat, nmax = 2)) {
    if (freq) weights[treat == i] <- weights[treat == i]/mean(weights[treat == i])
    else weights[treat == i] <- weights[treat == i]/sum(weights[treat == i])
  }

  ylab <- if (freq) "Count" else "Proportion"

  par(mfrow = c(2,2))
  # breaks <- pretty(na.omit(pscore), 10)
  breaks <- seq(minp, maxp, length = 11)
  xlim <- range(breaks)

  for (n in c("Raw Treated", "Matched Treated", "Raw Control", "Matched Control")) {
    if (startsWith(n, "Raw")) weights <- rep(1, length(treat))
    else weights <- x$weights

    if (endsWith(n, "Treated")) t <- 1
    else t <- 0

    if (freq) weights[treat == t] <- weights[treat == t]/mean(weights[treat == t])
    else weights[treat == t] <- weights[treat == t]/sum(weights[treat == t])

    #Create histogram using weights
    #Manually assign density, which is used as height of the bars. The scaling
    #of the weights above determine whether they are "counts" or "proportions".
    #Regardless, set freq = FALSE in plot() to ensure density is used for bar
    #height rather than count.
    pm <- hist(pscore[treat==t], plot = FALSE, breaks = breaks)
    pm[["density"]] <- vapply(seq_len(length(pm$breaks) - 1), function(i) {
      sum(weights[treat == t & pscore >= pm$breaks[i] & pscore < pm$breaks[i+1]])
    }, numeric(1L))
    plot(pm, xlim = xlim, xlab = xlab, main = n, ylab = ylab,
         freq = FALSE, col = "lightgray", ...)

    if (!startsWith(n, "Raw") && !is.null(q.cut)) abline(v = q.cut, lty=2)
  }

}
