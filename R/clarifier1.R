# AR1 test

phis = seq(0.99, 0.999, ((0.999 - 0.99)/5) )

N = 2^17
n = length(phis)
blocking = rep(0,n)
exact = rep(0,n)
ind = rep(0,n)
As = matrix(0,n,N)
Results = matrix(0,n,log2(N) )
ACFs = matrix(0,n,N+1)
for (k in 1:n) {
  print(k)
  As[k,] = arima.sim(model=list(order=c(1,0,0), ar=phis[k] ),n=N, rand.gen=rnorm,sd=1,mean=0)
  sigma2 = 1/( 1 - phis[k]^2 )
  ACFs[k,] = as.numeric(sigma2*ARMAacf(ar=c(phis[k]), ma=0, N));
  exact[k] = sigma2
  exact[k] = exact[k] + 2*sum( (1 - 1:( N-1 )/N)*ACFs[k,2:N ] )
  exact[k] = exact[k]/N
  
  # blocking
  x = As[k,]
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
  Results[k,] = ( s + cumsum( 2^(js-1)*gamma[js] ) )/(m-1)
  ind[k] = which.min( abs(exact[k] - Results[k,] ) )
  blocking[k] = Results[k,ind[k]]
}

blocking.eps = abs(blocking - exact)/exact

#plot(phis, blocking.eps )


# explorer code

K = 6
x = As[K,]
# num blockings
ind[K] - 1
m = length(x)
d = log2(length(x))
#ACF = c(acf(x,2^10, type="covariance")$acf)
Is = 0:(d-1)
Ns = 2^( d- Is )
Gamma = matrix(0, d, m +1 )
Gamma[1,] = ACFs[K,]
Sum = rep(0,d)
Sum[1 ] = sum(Gamma[1,])

for (k in 2:d ) {
  Gamma[k,1] = 0.5*Gamma[k-1, 1] + 0.5*Gamma[k-1, 2]
  for (i in 1:(Ns[k] - 1) ) {
    Gamma[k, i+1] = 0.25*Gamma[k-1, 2*(i) ] + 0.5*Gamma[k-1, 2*(i) + 1 ] + 0.25*Gamma[k-1, 2*(i)+2]
  }
  Sum[k ] = sum(Gamma[k,])
}

#error = cumsum(2^Is/(Ns - 1)*(    -exact[K]*(Ns + 1)  +  2*Sum/Ns  ) )
error = (    -exact[K]*(Ns + 1)  +  2*Sum/Ns  )/(Ns-1)
plot(0:(d-1), error )
#plot(0:(d-1), abs(error)/(sum(2^Is*Gamma[,2])) )
plot(0:(d-1), abs(Results[K,] - exact[K])/exact[K] )


nk = length(x)
y = rep(0,nk/2)
for (j in 1:nk/2) {
  y[j] = 0.5*(x[j] + x[nk/2+j])
}
x = y
ACF = c(acf(x,2^14, type="covariance")$acf)
sum(ACF)
k = k + 1
