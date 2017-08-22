data = read.table("/home/marius/Dokumenter/master/blocking/validation/build-teletubbies-Desktop_Qt_5_7_0_GCC_64bit-Release/data_1501957507286.txt")
data = read.table("/home/marius/Dokumenter/master/blocking/validation/build-teletubbies-Desktop_Qt_5_7_0_GCC_64bit-Release/data_1502047901965.txt")
data = read.table("/home/marius/Dokumenter/master/blocking/validation/build-teletubbies-Desktop_Qt_5_7_0_GCC_64bit-Release/data_1502834807876.txt")
data = read.table("/home/marius/Dokumenter/master/blocking/validation/build-teletubbies-Desktop_Qt_5_7_0_GCC_64bit-Release/data_1502916102276.txt")
new = (((data$V2-data$V1)/data$V1 )^2)^(1/4)
old = (((data$V3-data$V1)/data$V1 )^2)^(1/4)

new[new > 4.1] = NA
old[old > 4.1] = NA

explainatory = 2^20/(1/-log(data$V4))  #ar1

#ar2
phi1 = data$V4; phi2 = data$V5; 
x = -phi1/(2*phi2); y = -(- phi1^2 - 4*phi2)^0.5/(2*phi2);
z = (-phi2)^(-0.5); theta = atan2(y,x); b = atan( (phi1/(z*(1 - phi2)) - cos(theta) )/sin(theta) )
a = 1/(cos(b))
explainatory = 2^20*(-log(-phi2))/(2*(1 + log(abs(a) )))
explainatory[explainatory < 2e2] = NA

# ---------


old.fit = glm(old ~ log(explainatory) + I(log(explainatory)^2), family = poisson(link = sqrt) )
new.fit = glm(new ~ log(explainatory) + I(log(explainatory)^2), family = poisson(link = sqrt) )

#old.fit = glm(old ~ log(explainatory), family = poisson(link = log) )
#new.fit = glm(new ~ log(explainatory), family = poisson(link = log) )
#old.fit = glm(old ~ log(explainatory), family = inverse.gaussian(link = log), start=c(1.4720,-0.2676) )
new.fit = glm(new ~ log(explainatory) + I(log(explainatory)^2) )
old.fit = glm(old ~ log(explainatory) + I(log(explainatory)^2) )
summary(old.fit)
summary(new.fit)

newdata = data.frame(explainatory = seq(min(explainatory, na.rm=T),max(explainatory, na.rm=T),max(explainatory, na.rm=T)/10000 ) )
new.pred = predict(new.fit, newdata=newdata,se.fit = T, type="response")
old.pred = predict(old.fit, newdata=newdata,se.fit = T, type="response")

plot(explainatory,   new^2,pch=8,cex=0.05,col="blue", xlab = "obs/tau", ylab = "rel.error", log="x", ylim=c(0,3), xlim=c(200,1e5) )
points(explainatory, old^2,pch=8,cex=0.01,col="red" )
lines(newdata$explainatory, new.pred$fit^2, cex=1, lwd=1.5, col="blue")
#lines(newdata$explainatory, (new.pred$fit + 1.96*new.pred$se.fit)^2, lwd=2, lty=3, col="blue")
#lines(newdata$explainatory, (new.pred$fit - 1.96*new.pred$se.fit)^2, lwd=2, lty=3, col="blue")

lines(newdata$explainatory, old.pred$fit^2, cex=1, lwd=1.5, col="red")
#lines(newdata$explainatory, (pred2o$fit + 1.96*pred$se.fit), lwd=2, lty=3)
#lines(newdata$explainatory, (pred2o$fit - 1.96*pred$se.fit), lwd=2, lty=3)