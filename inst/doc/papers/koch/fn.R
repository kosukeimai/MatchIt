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
  bhat <- ols1$coefficients 
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
  legend(0.7, max(c(dx1$y, dx0$y)), lty = 1:2, col = 1:2, 
         legend = c("Treatment", "Control"), cex=0.5)
}

#model adjustment
ols.ate <- function(ols.m,sims=1000){
  b.sim <- as.matrix(draw.ols(ols.m,sims=sims))
  b.sim <- b.sim[,1:(ncol(b.sim)-1)]
  XX <- as.matrix(cbind(1,match.data(m1,"treat")[,m1$covariates]))
  y.hat <- XX%*%t(b.sim)
  ate <- matrix(0,sims,1)
  y.obs <- model.frame(ols.m)[,1]
  nn <- length(y.obs)
  if(nn==nrow(y.hat)){
    yy <- matrix(y.obs,nn,sims)
  }
  obj <- apply(yy-y.hat,2,mean)
  obj
}

#wrapper for zelig
zsim <- function(fml,dta){
  tt <- dta$dviswom
  z.out <- zelig(fml, data = dta, model="ls")
  x1 <- setx(z.out, data = dta[tt==1,],cond=T)
  x0 <- setx(z.out, data = dta[tt==0,],cond=T)
  x1[,dimnames(x.out)[[2]]=="dviswom"] <- 0
  x0[,dimnames(x.out)[[2]]=="dviswom"] <- 1
  s1 <- sim(z.out, x = x1, fn=NULL)
  s0 <- sim(z.out, x = x0, fn=NULL)
  ate <- c(s1$qi$ate.ev,-s0$qi$ate.ev)
  ate
}
