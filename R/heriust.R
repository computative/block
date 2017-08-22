n = 2^13
Sigma = matrix(0,n,n)
for (i in 1:n) {
  for (j in 1:n ) {
    t = abs(i-j)
    if(t == 0) {
      Sigma[i,i] = 1.2
    } else {
      Sigma[i,j] = exp(-t/2e2)*cos(t/3e1)
    }
  }
}
sdev = matrix(c(6,3,2,2),2,2)
X = rnorm(n,3,sdev)
for (k in 1:100) {
  X = rnorm(n,3,Sigma)
  a = acf(X, n, type="covariance", plot=FALSE)
  S = abs(sum(a$acf[a$acf < 0]))
  if(S < 50) {
    A = X
    print(S)
    acf(A, 400, type="covariance")
    break
  }
}

boot(A,statistic,100)
x = A

S = rep(0,12)1:040159
sdev = rep(0,12)
for (i in 1:12) {
  y = rep(0,length(x)/2)
  for (j in 1:length(x)/2) {
    y[j] = 0.5*(x[2*j-1] + x[2*j])
  }
  x = y
  N = length(x)
  c0 = var(x)*(N-1)/N
  S[i] = c0/(N-1)
  sdev[i] = sqrt(2/(N-1))*c0/(N-1)
}
plot(S, type="l", ylim=c(0,1e-8))
arrows(1:17, S-1.959964*sdev, 1:17, S+1.959964*sdev, length=0.01, angle=90, code=3, col="red")
