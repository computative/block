mydata = read.table("/home/marius/Dokumenter/fys4411/postmann.pat/code/dataget/data/data.txt") 
A = mydata$V1

x = A
d = log2(length(x))
m = 2^d
s = var(x)*(m-1)/m

for (i in 0:( d - 1) ) {
  nk = length(x)
  gamma = acf(x, 1, type="covariance", plot=FALSE)$acf[2]*nk/(nk-1)
  rho = acf(x, 1, type="correlation", plot=FALSE)$acf[2]*nk/(nk-1)
  if(abs(rho)<2/nk^0.5 ) {
    next
  }
  s = s + 2^i*gamma
  y = rep(0,nk/2)
  for (j in 1:nk/2) {
    y[j] = 0.5*(x[2*j-1] + x[2*j])
  }
  x = y
}

print(s/(m-1))


stat <- function(dat)
{
  mean(dat)
}

var(tsboot(A,stat, 2^8, l = 2^14, sim="geom")$t)
