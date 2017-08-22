data = read.table("../resources/data.txt") 

x = data$V1
nk = length(x)
d = log2(nk)
mu = mean(x)
s = matrix(0,d)
gamma = rep(0,d)
for (i in 0:(d-1) ) {
  nk = length(x)
  X = x - mu
  gamma[i+1] = sum( X[1:(nk-1)]*X[2:nk] )/(nk-1)
  s[i+1] = sum( X^2 )/nk
  x = 0.5*(x[2*(1:(nk/2))-1] + x[2*(1:(nk/2))])
}
s[d] = s[d-1]
ind = which.max(as.integer( rev( cumsum( rev( (gamma/s)^2*2^(d-0:(d-1) ) ) ) ) < qchisq(0.99, d:1) ) )
s[ind]/2^(d-(ind-1))
