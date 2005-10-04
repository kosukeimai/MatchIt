## Function for QQ summary stats
qqsum <- function (x, y){
  sx <- sort(x)
  sy <- sort(y)
  lenx <- length(sx)
  leny <- length(sy)
  if (leny < lenx) 
    sx <- approx(1:lenx, sx, n = leny)$y
  if (leny > lenx) 
    sy <- approx(1:leny, sy, n = lenx)$y
  dxy <- abs(sx-sy)
  meandiff <- mean(dxy)
  meddiff <- as.numeric(median(dxy))
  maxdiff <- max(dxy)
  invisible(list(meandiff=meandiff,
                 meddiff = meddiff,
                 maxdiff = maxdiff))
}
