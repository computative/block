stat <- function(dat)
{
  mean(dat)
}

boot.rho <- function(dat, inds)
{
  acf(x[inds],1, plot=F)$acf[2]
}


gamma = function(h,phi) 
{
  phi^abs(h)/(1-phi^2)
}

d = 20
k = 11
G = rep(0,2^(d-k) )
phi = 0.99
for (i in 1:2^(d-k) ) {
  print(i)
  for (j in 1:(2^(k+1)-1) ) {
    G[i] = G[i] + (2^k - abs(j - 2^k) )*gamma(2^k*(i-1) + j,phi)
  }
  G[i] = G[i]*2^(-2*k)
}
plot(G)
G

mydata = read.table("/home/marius/Dokumenter/fys4411/postmann.pat/code/dataget/data/data.txt") 
A = mydata$V1
A = A[1:2^14]
acf(A,15000)
#library(boot)
#boot(mydata$V1,statistic,500)
#boot(A,statistic,500)

var(tsboot(A, stat, 2^8, sim="geom", l=2^14)$t)[1]

# estimation 2
#A = arima.sim(model=list(ar=0.99),n=2^10)
x = A[1:2^14]
gamma = acf(x, length(x), type="covariance", plot=FALSE)$acf
S = 0
m = length(x)
for (t in 0:(m-1) ) {
  S = S + (m-t)*gamma[t+1]
}
print( S/(m^2) )


arima.sim(model=list(order=c(1,0,0), ar=0.9 ),n=2^14, rand.gen=rnorm,sd=1,mean=0)
acf(arima.sim(model=list(order=c(1,0,0), ar=0.9 ),n=2^14, rand.gen=rnorm,sd=1,mean=0),type="covariance")
lines(y,0.9^y/(1-0.9^2))

plot(arima.sim(model=list(order=c(1,0,0), ar=0.9 ),n=2^14, rand.gen=rnorm,sd=1,mean=3))
var(tsboot(A, stat, 2^9, sim="geom", l=10000)$t)[1]

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
    print(i-1)
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
print ( S/(m-1) )




if(false) {
  
  A = arima.sim(model=list(ar=0.999),n=2^16)
  x = A
  d = log2(length(x))
  m = 2^d
  S = var(x)*(m -1)/m
  s = 0
  gamma = acf(x, 2^14, type="covariance", plot=FALSE)$acf
  for (i in 0:( d -1) ) {
    gamma1 = 0
    for (j in 1:(2^(i+1) -1) ) {
      if (j+1 < 2^14) {
        gamma1 = gamma1 + 2^(-i)*gamma[j+1]*(2^i - abs(2^i - j))
      }
    }
    s = s + gamma1
  }
  print( (S+s)/(m-1) )
  
  
  
  A = arima.sim(model=list(ar=0.9),n=2^20)
  x = A
  d = log2(length(x))
  m = 2^d
  S = var(x)*(m -1)/m
  s = 0
  for (i in 0:( d -1 - 8) ) {
    nk = length(x)
    gamma1 = acf(x, 1, type="covariance", plot=FALSE)$acf[2]
    s = s + gamma1*(2^i)
    
    y = rep(0,nk/2)
    for (j in 1:nk/2) {
      y[j] = 0.5*(x[2*j-1] + x[2*j])
    }
    x = y
  }
  print( (S+s)/(m-1) )
}

nk = length(x)
gamma1 = acf(x/sqrt(m-1), 100, type="covariance")
abs(gamma1$acf[2])
y = rep(0,nk/2)
for (j in 1:nk/2) {
  y[j] = 0.5*(x[2*j-1] + x[2*j])
}
x = y

