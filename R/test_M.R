statboot <- function(X,inds)
{
  dat = X[inds]
  S = 0
  for (i in k:( d - 1) ) {
    nk = length(dat)
    gamma1 = acf(dat, 1, type="covariance", plot=FALSE)$acf[2]#*nk/(nk-1)
    S = S + 2^i*gamma1
    
    y = rep(0,nk/2)
    for (j in 1:nk/2) {
      y[j] = 0.5*(dat[2*j-1] + dat[2*j])
    }
    dat = y
  }
  S
}

statacf <- function(X,inds)
{
  x = X[inds]
  c(acf(x,1,type="correlation",plot=F)$acf[2] )
}

statci1 <- function(X,inds)
{
  x = X[inds]
  n = length(x)
  sum( x[1:(n-1)]*x[2:n] )/(n-1) - mean(x)^2
}


statacfcor <- function(X,inds)
{
  nk = length(X)
  acf(X[inds], 1, type="covariance", plot=FALSE)$acf[2]*nk/(nk-1)
}

rhoboot <- function(X,inds)
{
  dat = X[inds]
  S = 0
  for (i in k:( d - 1) ) {
    nk = length(dat)
    rho1 = acf(dat, 1, plot=FALSE)$acf[2]*nk/(nk-1)
    S = S + rho1
    y = rep(0,nk/2)
    for (j in 1:nk/2) {
      y[j] = 0.5*(dat[2*j-1] + dat[2*j])
    }
    dat = y
  }
  S
}

ts.statmean <- function(X)
{
  mean(X)
}

statmean <- function(X,inds)
{
  mean(X[inds])
}

library(boot)



x = A
acf(x,1000, type="covariance")
m = length(x)
d = log2(length(x))
k = 0

  nk = length(x)
  x = 0.5*(x[2*(1:(nk/2))-1] + x[2*(1:(nk/2))])
  acf(x,100, type="correlation")
  k = k + 1
  

stat <- function(dat, inds)
{
  S = 0
  x = dat[inds]
  for (i in k:( d - 1) ) {
    nk = length(x)
    rho = acf(x, 1, type="correlation", plot=FALSE)$acf[2]
    S = S + rho^2*nk
    x = 0.5*(x[2*(1:(nk/2))-1] + x[2*(1:(nk/2))])
  }
  S
}

#star = c(tsboot(z, stat, 2^9, sim="geom", l=length(z) )$t)
z = x
star = c(boot(z, stat, 2^12)$t)
star = c(boot(z, statci1, 2^14)$t)
mean(star)
sd(star)
star = c(boot(z, statacf, 2^14 )$t)
var(star*sqrt(length(z)),na.rm=T)
mean(star,na.rm=T)
sd(star)
Is = 0:(d-1-k)
sum(  2^(Is+k)*(2^(k+1-d)*var(z) - varmu*(2^(d-Is-k) +1 ) )/(2^(d-Is-k) - 1) )
var(star)
sd(star)
hist( star, 50, prob=T)
ts = seq(-1500,1500,5)
fts = dnorm(ts,0, (s1)^0.5 ) 
lines(ts,fts)