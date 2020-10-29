jitter.pscore <- function(x, interactive, pch = 1, ...){

  .pardefault <- par(no.readonly = TRUE)
  on.exit(par(.pardefault))

  treat <- x$treat
  pscore <- x$distance
  s.weights <- if (is.null(x$s.weights)) rep(1, length(treat)) else x$s.weights
  weights <- x$weights * s.weights
  matched <- weights > 0
  q.cut <- x$q.cut
  jitp <- jitter(rep(1,length(treat)), factor=6)+(treat==1)*(weights==0)-(treat==0) - (weights==0)*(treat==0)
  cswt <- sqrt(s.weights)
  cwt <- sqrt(weights)
  minp <- min(pscore, na.rm = TRUE)
  maxp <- max(pscore, na.rm = TRUE)

  plot(pscore, xlim = c(minp - 0.05*(maxp-minp), maxp + 0.05*(maxp-minp)), ylim = c(-1.5,2.5),
       type="n", ylab="", xlab="Propensity Score",
       axes=F,main="Distribution of Propensity Scores",...)
  if (!is.null(q.cut)) abline(v = q.cut, col = "grey", lty = 1)

  #Matched treated
    points(pscore[treat==1 & matched], jitp[treat==1 & matched],
           pch = pch, cex = cwt[treat==1 & matched], ...)
    #Matched control
    points(pscore[treat==0 & matched], jitp[treat==0 & matched],
           pch = pch, cex = cwt[treat==0 & matched], ...)
    #Unmatched treated
    points(pscore[treat==1 & !matched], jitp[treat==1 & !matched],
           pch = pch, cex = cswt[treat==1 & matched],...)
    #Unmatched control
    points(pscore[treat==0 & !matched], jitp[treat==0 & !matched],
           pch = pch, cex = cswt[treat==0 & matched], ...)

  axis(1)

  center <- mean(par("usr")[1:2])
  text(center, 2.5, "Unmatched Treated Units", adj = .5)
  text(center, 1.5, "Matched Treated Units", adj = .5)
  text(center, 0.5, "Matched Control Units", adj = .5)
  text(center, -0.5, "Unmatched Control Units", adj = .5)
  box()

  if (interactive) {
    print("To identify the units, use first mouse button; to stop, use second.")
    identify(pscore, jitp, names(treat), atpen = TRUE)
  }
}
