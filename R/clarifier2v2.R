library(boot)

# AR2 test

n = 6080
N = 2^20
post = floor(log10(N))
pre = 10^(log10(N) - post)
phi2s = -exp(-10^(runif( n,min=2,max=6 ) - post )/pre )
phi1s = rep(0,n)
for (i in 1:n) {
  phi1s[i] = runif(1, 0,(-4*phi2s[i] )^0.5 )
}

Blocking = rep(0,n)
blocking = rep(0,n)
exact = rep(0,n)
strap = rep(0,n)
ind = rep(0,n)
As = rep(0,N)
Results = matrix(0,n,log2(N) )
ACFs = rep(0,N)
s = matrix(0,n,log2(N))
gamma = matrix(0,n,d)
for (k in 1:n) {
  print(k)
  As = arima.sim(list(ar=c(phi1s[k],phi2s[k])), n = N, rand.gen=rgamma,shape=1,scale=1)
  sigma2 = (1 - phi2s[k])/( (1 + phi2s[k])*( (1-phi2s[k])^2 - phi1s[k]^2 ) )
  ACFs = as.numeric(sigma2*ARMAacf(ar=c(phi1s[k],phi2s[k]), ma=0, N-1));
  exact[k] = sigma2
  exact[k] = exact[k] + 2*sum( (1 - 1:( N-1 )/N)*ACFs[2:N ] )
  exact[k] = exact[k]/N
  
  # bootstrapping
  strap[k] = var(tsboot(As,mean,R=8, sim="geom", l=N)$t)
  
  # blocking
  #x = As
  #nk = length(x)
  #d = log2(nk)
  #mu = mean(x)
  #for (i in 0:(d-1) ) {
  #  nk = length(x)
  #  X = x - mu
  #  gamma[k,i+1] = sum( X[1:(nk-1)]*X[2:nk] )/(nk-1)
  #  s[k,i+1] = sum( X^2 )/nk
  #  x = 0.5*(x[2*(1:(nk/2))-1] + x[2*(1:(nk/2))])
  #}
  #s[k,d] = s[k,d-1]
  #ind[k] = which.max(as.integer( rev( cumsum( rev( (gamma[k,]/s[k,])^2*2^(d-0:(d-1) ) ) ) ) < qchisq(0.99, d:1) ) )
  #blocking[k] =  s[k,ind[k]]/2^(d-(ind[k]-1)) + sum(0.85^((ind[k]:d - 1) )*gamma[k,ind[k]:d])/2^(d-(ind[k]-1))
  #Blocking[k] =  s[k,ind[k]]/2^(d-(ind[k]-1))
}

#Blocking.eps = ((Blocking - exact)/exact)^2
#blocking.eps = ((blocking - exact)/exact)^2
strap.eps = ((strap-exact)/exact)^2

x = -phi1s/(2*phi2s); y = -(- phi1s^2 - 4*phi2s)^0.5/(2*phi2s);
z = (-phi2s)^(-0.5); theta = atan2(y,x); b = atan( (phi1s/(z*(1 - phi2s)) - cos(theta) )/sin(theta) )
a = 1/(cos(b))
expl = N*(-log(-phi2s))/(2*(1 + log(abs(a) )))

plot(expl, strap.eps, log = "xy" ,pch=8,cex=0.5,xlim=c(1e1,1e6),ylim=c(0.001,2) )
points(expl, blocking.eps, pch=1 ,cex=.5)
points(expl, Blocking.eps, pch=3 ,cex=.5)

# explorer code

K = 2
x = As[K,]
# num blockings
ind[K] - 1
#ACF = c(acf(x,2^14, type="covariance")$acf)
sum(ACF)
m = length(x)
d = log2(length(x))

Is = 0:(d-1)
Ns = 2^( d- Is )
Gamma = matrix(0, d, m )
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
terms = (    -exact[K]*(Ns + 1)  +  2*Sum/Ns  )/(Ns-1)
error = cumsum( 2^Is*(    -exact[K]*(Ns + 1)  +  2*Sum/Ns  )/(Ns-1) )
rest = rev(cumsum(rev(Gamma[,2])*2^rev(Is) ) )

plot(0:(d-1), abs((error ) /(m-1))/exact[K] )
plot(0:(d-1), abs((rest ) /(m-1))/exact[K] )
plot(0:(d-1), abs((error - rest ) /(m-1))/exact[K] )

plot(0:(d-1), abs(terms) )
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



