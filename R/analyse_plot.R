options(latexcmd="/usr/local/texlive/2017/bin/x86_64-linux/pdflatex")
options(latexcmd="latex")
options(tikzLatex="pdflatex")
library(tikzDevice)





tikz("/home/marius/Dokumenter/master/blocking/bilder/fig1.tex", width = 6.5, height = 3.5, pointsize=11)
par(mfrow=c(1,2))
par(mgp=c(2,0.5,0), mar=c(5,4,4,2), cex=1.24)
par(ps=8, cex=1, cex.main=1.3, cex.lab=1.3, cex.axis=1.3, cex.sub=1, font.main=1)
xs = seq(2,6,0.05)

data = read.table("/home/marius/Dokumenter/master/blocking/validation/build-teletubbies-Desktop_Qt_5_7_0_GCC_64bit-Release/AR1Log2n24_56245.txt"); n = 2^24 # G 0.99
old = (((data$V3-data$V1)/data$V1 )^2)
phi = data$V4; explainatory = n/(1/-log(phi))  #ar1

explainatory = explainatory + runif(length(explainatory), 0,20)
explainatory[explainatory < 1e2] = NA

par(mar=c(3.6,3.6,2,0.1))

old.fit = glm(old ~ log(explainatory), family = Gamma(link = log) )
summary(old.fit)
newdata = data.frame(explainatory = 10^xs )
old.pred = predict(old.fit, newdata=newdata,se.fit = T, type="response")
plot(explainatory ,   old,pch="*",cex=0.3, ylab="Relative error square",
      xlab = "$n/\\tau$", log = "x",xlim=c( 1e2,1e6) ,ylim=c(0, 0.5),main="Error AR(1) Gamma process" )
lines(newdata$explainatory, old.pred$fit, cex=1, lwd=1.5, col="black")
lines(newdata$explainatory, (old.pred$fit + 1.96*old.pred$se.fit), lwd=2, lty=3)
lines(newdata$explainatory, (old.pred$fit - 1.96*old.pred$se.fit), lwd=2, lty=3)


par(mar=c(3.6,2,2,1.7))
data = read.table("/home/marius/Dokumenter/master/blocking/validation/build-teletubbies-Desktop_Qt_5_7_0_GCC_64bit-Release/AR2Log2n24_46310.txt"); n = 2^24 # G 0.99
old = (((data$V3-data$V1)/data$V1 )^2)
phi1 = data$V4; phi2 = data$V5; 
x = -phi1/(2*phi2); y = -(- phi1^2 - 4*phi2)^0.5/(2*phi2);
z = (-phi2)^(-0.5); theta = atan2(y,x); b = atan( (phi1/(z*(1 - phi2)) - cos(theta) )/sin(theta) )
a = 1/(cos(b))
explainatory = n*(-log(-phi2))/(2*(1 + log(abs(a) )))
explainatory = explainatory + runif(length(explainatory), 0,20)
explainatory[explainatory < 1e2] = NA

old.fit = glm(old ~ log(explainatory), family = Gamma(link = log) )
summary(old.fit)
newdata = data.frame(explainatory = 10^xs )
old.pred = predict(old.fit, newdata=newdata,se.fit = T, type="response")
plot( (explainatory + runif(length(explainatory), 0,20)),   old,pch="*",cex=0.3,xlab = "$n/\\tau$", 
      log = "x",xlim=c( 1e2,1e6) ,ylim=c(0, 1) , ylab=NA,main="Error AR(2) Normal process")
lines(newdata$explainatory, old.pred$fit, cex=1, lwd=1.5, col="black")
lines(newdata$explainatory, (old.pred$fit + 1.96*old.pred$se.fit), lwd=2, lty=3)
lines(newdata$explainatory, (old.pred$fit - 1.96*old.pred$se.fit), lwd=2, lty=3)

dev.off()
