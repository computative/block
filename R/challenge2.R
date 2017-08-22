# AR1 test

stat <- function(dat)
{
  mean(dat)
}

gamma = function(h,phi) 
{
  phi^abs(h)/(1-phi^2)
}

library(boot)

n = 5
phi2s = -seq(0.75, 0.999, ((0.999 - 0.75)/(n-1)) )
phi1s = rep(0,n)

for (i in 1:n) {
  phi1s[i] = runif(1, 0,(-4*phi2s[i] )^0.5 )
}

N = 2^20
blocking = rep(0,n)
bootstrap = rep(0,n)
exact = rep(0,n)

for (k in 1:n) {
  print(k)
  A = arima.sim(list(ar=c(phi1s[k],phi2s[k])), n = N, rand.gen=rnorm,sd=1,mean=0)
  sigma2 = var(A)
  ACF = sigma2*ARMAacf(ar=c(phi1s[k],phi2s[k]), ma=0, N);
  exact[k] = sigma2
  exact[k] = exact[k] + 2*sum( (1 - 1:( N-1 )/N)*ACF[2:N ] )
  exact[k] = exact[k]/N
  
  z = c(-1, phi1s[k],phi2s[k])
  a = abs(polyroot(z)[1])
  tau = 2^4*max(100, ( 1 - log(a) )/log(a) )
  tau = min(tau,N)
  
  # bootstrap
  bootstrap[k] = var(tsboot(A, stat, 2^8, sim="geom", l=max(abs(tau),50) )$t)[1]
  
  # blocking
  x = A
  d = log2(length(x))
  m = 2^d
  s = var(x)*(m-1)/m
  gamma = rep(0,d)
  rho = rep(0,d)
  for (i in 0:( d - 1) ) {
    nk = length(x)
    temp = acf(x, 1, type="covariance", plot=FALSE)$acf
    gamma[i+1] = temp[2]*nk/(nk-1);
    rho[i+1] = temp[2]/temp[1]*nk/(nk-1); 
    y = rep(0,nk/2)
    for (j in 1:nk/2) {
      y[j] = 0.5*(x[2*j-1] + x[2*j])
    }
    x = y
  }
  
  
  observable = rev(cumsum(rev(rho)))
  iterate = abs(observable + 0.725) > 1.34/8
  
  j = 1
  while(iterate[j]) {
    s = s + gamma[j+1]
    j = j + 1
  }
  blocking[k] = s/(m-1)
}

blocking.eps = abs(blocking - exact)/exact
bootstrap.eps = abs(bootstrap - exact)/exact

plot(phi1s, bootstrap.eps, ylim = c(0.001, 20 ) )#, log="y" )
points(phi1s, blocking.eps ,pch=2)
