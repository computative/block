mydata = read.table("/home/marius/Dokumenter/fys4411/postmann.pat/code/dataget/data/data.txt") 

x = mydata$V1
n = length(x)
d = log2(n)
s = rep(0,d)
gamma = rep(0,d)
mu = mean(x)
for (i in 0:(d-1)) {
  nk = length(x)
  X = x - mu
  gamma[i+1] = sum( X[1:(nk-1)]*X[2:nk] )/(nk-1)
  s[i+1] = var(x) #sum( X^2 )/nk
  x = 0.5*(x[2*(1:(nk/2))-1] + x[2*(1:(nk/2))])
}
s[d] = s[d-1]
ind = which.max(as.integer(gamma/s < 1.96/sqrt(2^(d-0:19))))
#ind = which.max(as.integer( rev( cumsum( rev( (gamma/s)^2*2^(d-0:19) ) ) ) < qchisq(0.95, 20:1) ) )
print(s[ind]/2^(d-(ind-1) ) )