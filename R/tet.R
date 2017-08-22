# AR1 test

phi = 0.999
N = 2^10
A = arima.sim(model=list(order=c(1,0,0), ar=phi ),n=N, rand.gen=rnorm,sd=1,mean=0)
n = length(A)
sigma2 = 1/( 1 - phi^2 )
ACF = sigma2*as.numeric(ARMAacf(ar=c(phi), ma=0, N)) ;
varX = sigma2
varX = varX + 2*sum( (1 - 1:( N-1 )/N)*ACF[2:N ] )
varX = varX/N

E.diff = varX*(1 - 2*n/(n-1) ) + 2/(n*(n-1))*sum(ACF)
print(E.diff)
s = 0
for (i in 1:2000) {
  A = arima.sim(model=list(order=c(1,0,0), ar=phi ),n=N, rand.gen=rnorm,sd=1,mean=0)
  s = s + acf(A,1,type="covariance")$acf[2]*n/(n-1) - ACF[2]
}
print(s/2000)










############

# AR1 test

n = 5
phi2s = -seq(0.75, 0.999, ((0.999 - 0.75)/(n-1)) )
phi1s = rep(0,n)
d = 10

for (i in 1:n) {
  phi1s[i] = runif(1, 0,(-4*phi2s[i] )^0.5 )
}

N = 2^d
S = 0
K = 3
n = 10000
for (o in 1:n) {
  x = arima.sim(list(ar=c(phi1s[K],phi2s[K])), n = N, rand.gen=rnorm,sd=1,mean=0)
  mu = mean(x)
  S = S - 1/n*sum( (x - mu)^2  )
  for (i in 0:(d-1))
  {
    nk = length(x)
    js = 1:(nk/2)
    gamma1 = 2/n*sum( (x[js] - mu)*(x[js+nk/2] - mu)  )
    S = S + 2^i*gamma1
    y = rep(0,nk/2)
    for (j in 1:nk/2) {
      y[j] = 0.5*(x[2*j-1] + x[2*j])
    }
    x = y
  }
}
print(S/n)
