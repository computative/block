mydata = read.table("/home/marius/Dokumenter/fys4411/postmann.pat/code/dataget/data/data.txt") 

x = mydata$V1
n = length(x)
d = log2(n)
s = var(x)*(n-1)/n
mu = mean(x)
J = d-1
for (i in 0:(d-1) ) {
  nk = length(x)
  X = x - mu
  #gamma = sum( x[1:(nk-1)]*x[2:nk] )/(nk-1) - mu^2
  gamma = sum( X[1:(nk-1)]*X[2:nk] )/(nk-1)
  Vx = sum( X^2 )/nk
  x = 0.5*(x[2*(1:(nk/2))-1] + x[2*(1:(nk/2))])
  if (abs(gamma) < 1.96*Vx/(nk)^0.5) {
    J = i
    break
  }
  s = s + 2^i*gamma
}
if (J == (d-1) ) {
  s = s/(n-1)
} else {
  s = s/(n-2^(J+1))
}

print(s)