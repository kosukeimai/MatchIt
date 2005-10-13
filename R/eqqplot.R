eqqplot <- function(x, y, plot.it = TRUE, xlab = deparse(substitute(x)),
                    ylab = deparse(substitute(y)), ...)
{ ## empirical quantile-quantile plot; hacked from qqplot() in stats.
  sx <- sort(x)
  sy <- sort(y)
  lenx <- length(sx)
  leny <- length(sy)
  if (leny < lenx)
    sx <- approx(1:lenx, sx, n = leny, method = "constant")$y
  if (leny > lenx)
    sy <- approx(1:leny, sy, n = lenx, method = "constant")$y
  if (plot.it)
    plot(sx, sy, xlab = xlab, ylab = ylab, ...)
  invisible(list(x = sx, y = sy))
}
