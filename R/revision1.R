# AR2 test

n = 5
phis = seq(0.99, 0.9999, ((0.9999 - 0.99)/5) )

N = 2^20
blocking = rep(0,n)
oblocking = rep(0,n)
exact = rep(0,n)
ACFs = matrix(0,n,N+1)
for (k in 1:n) {
  print(k)
  A = arima.sim(model=list(order=c(1,0,0), ar=phis[k] ),n=N, rand.gen=rnorm,sd=1,mean=0)
  sigma2 = 1/( 1 - phis[k]^2 )
  ACFs[k,] = as.numeric(sigma2*ARMAacf(ar=c(phis[k]), ma=0, N));
  exact[k] = sigma2
  exact[k] = exact[k] + 2*sum( (1 - 1:( N-1 )/N)*ACFs[k,2:N ] )
  exact[k] = exact[k]/N
  
  # blocking
  x = A
  d = log2(length(x))
  m = 2^d
  s = var(x)*(m-1)/m
  S = s
  mu = mean(x)
  for (i in 0:( d - 1) ) {
    nk = length(x)
    X = x - mu
    gamma = sum( X[1:(nk-1)]*X[2:nk] )/(nk-1)
    Vx = sum( X^2 )/nk
    oblocking[k] = var(x)/nk
    x = 0.5*(x[2*(1:(nk/2))-1] + x[2*(1:(nk/2))])
    if (abs(gamma) < 1.96*Vx/(nk)^0.5) {
      next
    }
    s = s + 2^i*gamma
  }
  blocking[k] = abs(s)/(N-1)
}
oblocking.eps = abs(oblocking - exact)/exact
blocking.eps = abs(blocking - exact)/exact
blocking
oblocking
exact
plot(oblocking.eps, pch=2, ylim=c(0,max(oblocking.eps,blocking.eps) ))
points(blocking.eps)