#simulating posterior draws
draw.ols <- function(ols1,sims=1000){
  rinvchisq <- function (n, df, scale = 1/df){
    if ((length(scale) != 1) & (length(scale) != n))
      stop("scale should be a scalar or a vector of the same length as x")
    if (df <= 0)
      stop("df must be greater than zero")
    if (any(scale <= 0))
      stop("scale must be greater than zero")
    return((df * scale)/rchisq(n, df = df)) }
  s2 <- summary(ols1)$sigma^2 
  bhat <- na.omit(ols1$coefficients)
  k <- length(bhat)
  n <- nrow(model.matrix(ols1))
  vhat <- summary(ols1)$cov.unscaled 
  sim.sig2 <- rinvchisq(sims,n-k,s2)
  sim.beta <- matrix(0,sims,k)
  for(i in 1:sims){
    sim.beta[i,] <- rmvnorm(1,mean=bhat,sigma=vhat*sim.sig2[i])
  }
  sim.beta <- as.data.frame(cbind(sim.beta,sim.sig2))
  names(sim.beta) <- c(names(bhat),"sig2")
  return(sim.beta)
}
#overlaying densities
doverlay <- function(x1, x0, xlab = "", main = "", lines = FALSE, ...)
{
  x <- c(x1,x0)
  minobs <- min(x)
  maxobs <- max(x)
  dx1 <- density(x1, from = minobs, to = maxobs)
  dx0 <- density(x0, from = minobs, to = maxobs)
    matplot(dx0$x, cbind(dx1$y, dx0$y), type = "l", ylab = "Density", 
            xlab = xlab, main = main, ...)
  legend(-0.2, max(c(dx1$y, dx0$y)), lty = 1:2, col = 1:2, 
         legend = c("Full Data", "Matched Data"), cex=0.5)
}

#model adjustment
ols.ate <- function(ols.m,sims=1000,treat){
  b.sim <- as.matrix(draw.ols(ols.m,sims=sims))
  b.sim <- b.sim[,1:(ncol(b.sim)-1),drop=F]
  dta <- eval(ols.m$call$data,parent.frame())
  tt <- dta[,treat]
  dta[,treat] <- 1-dta[,treat]
  xx <- model.matrix(eval(ols.m$call$formula),data=dta)
  xx <- xx[,!is.na(coef(ols.m))]
  y.hat <- xx%*%t(b.sim)
  ate <- matrix(0,sims,1)
  y.obs <- model.frame(ols.m)[,1]
  nn <- length(y.obs)
  if(nn==nrow(y.hat)){
    yy <- matrix(y.obs,nn,sims)
  } else {
    stop("Dimensions of yobs and yhat off")
  }
  ss <- (1-tt)*y.hat+tt*yy-(tt*y.hat+(1-tt)*yy)
  obj <- apply(ss,2,mean)
  obj
}

#dta <- dta.full
#dta <- dta.match
#wrapper for zelig
zsim <- function(fml,dta,num=1000,zz=T,treat="dviswom"){
  if(zz){
    tt <- dta$dviswom
    z.out <- zelig(fml, data = dta, model="ls")
    dta$dviswom <- 1-dta$dviswom
    x1 <- setx(z.out, data = dta[tt==1,],cond=T)
    x0 <- setx(z.out, data = dta[tt==0,],cond=T)
    s1 <- sim(z.out, x = x1, fn=NULL, num=1000)
    s0 <- sim(z.out, x = x0, fn=NULL, num=1000)
    ate <- c(s1$qi$ate.ev,-s0$qi$ate.ev)   
  } else {
    ols.m <- lm(fml,data=dta)
    ate <- ols.ate(ols.m,sim=num,treat)
  }
  return(ate)
}
