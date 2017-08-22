# AR1 test

stat <- function(dat)
{
  mean(dat)
}

library(boot)

gamma = function(h,phi) 
{
  phi^abs(h)/(1-phi^2)
}

phi = 0.999
N = 2^18
A = arima.sim(model=list(order=c(1,0,0), ar=phi ),n=N, rand.gen=rnorm,sd=1,mean=0)
m = length(A)
exact = 0
for ( t in (1-m):(m-1) ) {
  exact = exact + (1 - t/m)*gamma(t,phi)
}
exact = exact/m
var(tsboot(A, stat, 2^8,l = 4500, sim="geom")$t)

x = A
d = log2(length(x))
S = var(x)*(m - 1)/m
mu = mean(x)
for (i in 0:( d - 1) ) {
  print(i)
  nk = length(x)
  js = 1:(nk/2)
  gamma1 = 1/(n-1)*sum( (x[js] - mu)*(x[js+1] - mu)  )
  if ( gamma1 < 0 ) {
    break
  }
  S = S + gamma1*(2^i)
  y = rep(0,nk/2)
  for (j in 1:nk/2) {
    y[j] = 0.5*(x[2*j-1] + x[2*j] )
  }
  x = y
}
print( S/(m-1) )






library(boot)
phis = seq(0.99, 0.999, ((0.999 - 0.99)/10) )

N = 2^15
n = length(phis)
blocking = rep(0,n)
bootstrap = rep(0,n)
exact = rep(0,n)
for (k in 1:n) {
  print(k)
  A = arima.sim(model=list(order=c(1,0,0), ar=phis[k] ),n=N, rand.gen=rnorm,sd=1,mean=0)
  
  # bootstrap
  tau = -2^4/log(phis[k])
  bootstrap[k] = var(tsboot(A, stat, 2^8, sim="geom", l=max(abs(tau),50) )$t)[1]
  
  # blocking
  x = A
  d = log2(length(x))
  m = 2^d
  S = var(x)*(m -1)/m
  gamma1 = S
  gamma1p = 0
  for (i in 0:( d - 1) ) {
    nk = length(x)
    gamma1p = gamma1
    gamma1 = acf(x, 1, type="covariance", plot=FALSE)$acf[2]*nk/(nk-1)
    if (abs(gamma1) > abs(gamma1p)) {
      print( S/(m-1) )
      print( (2^d - 2^(i-1) )*abs(gamma1p)/(m-1) )
      break
    }
    S = S + gamma1*(2^i)
    
    y = rep(0,nk/2)
    for (j in 1:nk/2) {
      y[j] = 0.5*(x[2*j-1] + x[2*j])
    }
    x = y
  }
  blocking[k] = S/(m-1)
  
  # exact
  exact[k] = 0
  for ( t in (1-m):(m-1) ) {
    exact[k] = exact[k] + (1 - t/m)*gamma(t,phis[k])
  }
}
exact = exact/m

blocking.eps = abs(blocking - exact)/exact
bootstrap.eps = abs(bootstrap - exact)/exact

plot(phis, bootstrap.eps, ylim = c(0.001, 2 ) )#, log="y" )
points(phis, blocking.eps ,pch=2)
