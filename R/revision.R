# AR2 test

n = 5
phi2s = -seq(0.75, 0.999, ((0.999 - 0.75)/(n-1)) )
phi1s = rep(0,n)

for (i in 1:n) {
  phi1s[i] = runif(1, 0,(-4*phi2s[i] )^0.5 )
}

N = 2^20
oblocking = rep(0,n)
blocking = rep(0,n)
exact = rep(0,n)
ACFs = matrix(0,n,N+1)
for (k in 1:n) {
  print(k)
  A = arima.sim(list(ar=c(phi1s[k],phi2s[k])), n = N, rand.gen=rnorm,sd=1,mean=0)
  sigma2 = (1 - phi2s[k])/( (1 + phi2s[k])*( (1-phi2s[k])^2 - phi1s[k]^2 ) )
  ACFs[k,] = as.numeric(sigma2*ARMAacf(ar=c(phi1s[k],phi2s[k]), ma=0, N));
  exact[k] = sigma2
  exact[k] = exact[k] + 2*sum( (1 - 1:( N-1 )/N)*ACFs[k,2:N ] )
  exact[k] = exact[k]/N
  
  # blocking
  x = A
  d = log2(length(x))
  m = 2^d
  s = var(x)*(m-1)/m
  mu = mean(x)
  penalty = 0
  J = d-1
  for (i in 0:( d - 1) ) {
    nk = length(x)
    X = x - mu
    gamma = sum( x[1:(nk-1)]*x[2:nk] )/(nk-1)- mu^2
    #gamma = sum( X[1:(nk-1)]*X[2:nk] )/(nk-1)
    Vx = sum( X^2 )/nk
    oblocking[k] = var(x)/nk
    x = 0.5*(x[2*(1:(nk/2))-1] + x[2*(1:(nk/2))])
    if (abs(gamma) < 1.96*Vx/(nk)^0.5) {
      #next
      penalty = penalty + 1
      if (penalty >= 2) {
        J = i
        break
      }
      #next
      #break
    }
    s = s + 2^i*gamma
  }
  if (J == (d-1) ) {
    blocking[k] = abs(s)/(N-1)
  } else {
    blocking[k] = abs(s)/(N-2^(J+1))
  }
}
oblocking.eps = abs(oblocking - exact)/exact
blocking.eps = abs(blocking - exact)/exact
oblocking
blocking
exact
plot(oblocking.eps, pch=2, ylim=c(0,max(oblocking.eps,blocking.eps) ) )
points(blocking.eps)
sum(oblocking.eps)
sum(blocking.eps)
