statboot <- function(X,inds)
{
  dat = X[inds]
  S = 0
  for (i in k:( d - 1) ) {
    nk = length(dat)
    gamma1 = acf(dat, 1, type="covariance", plot=FALSE)$acf[2]*nk/(nk-1)
    S = S + gamma1/var(y)
    
    y = rep(0,nk/2)
    for (j in 1:nk/2) {
      y[j] = 0.5*(dat[2*j-1] + dat[2*j])
    }
    dat = y
  }
  S
}

library(boot)



x = A
acf(x,1000, type="covariance")
m = length(x)
d = log2(length(x))
k = 0

  nk = length(x)

  y = rep(0,nk/2)
  for (j in 1:nk/2) {
    y[j] = 0.5*(x[2*j-1] + x[2*j])
  }
  x = y
  acf(x,100, type="covariance")
  k = k + 1
  

stat <- function(dat)
{
  S = 0
  for (i in k:( d - 1) ) {
    nk = length(dat)
    gamma1 = acf(dat, 1, type="covariance", plot=FALSE)$acf[2]*nk/(nk-1)
    S = S + gamma1/var(y)
    
    y = rep(0,nk/2)
    for (j in 1:nk/2) {
      y[j] = 0.5*(dat[2*j-1] + dat[2*j])
    }
    dat = y
  }
  S
}

#star = c(tsboot(z, stat, 2^9, sim="geom", l=length(z) )$t)
z = x
star = c(boot(z, statboot, 2^9 )$t)
mean(star)
hist(star+0.725,20,prob=T)
