# AR1 test

n = 5
phi2s = -seq(0.75, 0.999, ((0.999 - 0.75)/(n-1)) )
phi1s = rep(0,n)

for (i in 1:n) {
  phi1s[i] = runif(1, 0,(-4*phi2s[i] )^0.5 )
}

N = 2^18
blocking = rep(0,n)
exact = rep(0,n)

for (k in 1:n) {
  print(k)
  A = arima.sim(list(ar=c(phi1s[k],phi2s[k])), n = N, rand.gen=rnorm,sd=1,mean=0)
  sigma2 = (1 - phi2s[k])/( (1 + phi2s[k])*( (1-phi2s[k])^2 - phi1s[k]^2 ) )
  ACF = sigma2*ARMAacf(ar=c(phi1s[k],phi2s[k]), ma=0, N);
  exact[k] = sigma2
  exact[k] = exact[k] + 2*sum( (1 - 1:( N-1 )/N)*ACF[2:N ] )
  exact[k] = exact[k]/N
  
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
  js = 1:d
  results = ( s + cumsum( 2^(js-1)*gamma[js] ) )/(m-1)
  ind = which.min( abs(exact[k] -results ) )
  blocking[k] = results[ind]
}

blocking.eps = abs(blocking - exact)/exact

plot(phi1s, blocking.eps )
