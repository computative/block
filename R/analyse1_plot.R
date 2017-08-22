data = read.table("/home/marius/Dokumenter/master/blocking/validation/build-teletubbies-Desktop_Qt_5_7_0_GCC_64bit-Release/data_1503095604092.txt"); n = 2^20 # G AR1 0.99
data = read.table("/home/marius/Dokumenter/master/blocking/validation/build-teletubbies-Desktop_Qt_5_7_0_GCC_64bit-Release/data_1503095661261.txt"); n = 2^20 # G AR2 0.99
data = read.table("/home/marius/Dokumenter/master/blocking/validation/build-teletubbies-Desktop_Qt_5_7_0_GCC_64bit-Release/data_1503096488283.txt"); n = 2^20 # N AR1 0.99
data = read.table("/home/marius/Dokumenter/master/blocking/validation/build-teletubbies-Desktop_Qt_5_7_0_GCC_64bit-Release/data_1503096540774.txt"); n = 2^20 # N AR2 0.99
data = read.table("/home/marius/Dokumenter/master/blocking/validation/build-teletubbies-Desktop_Qt_5_7_0_GCC_64bit-Release/data_1503096833832.txt"); n = 2^24 # N AR1 0.99
data = read.table("/home/marius/Dokumenter/master/blocking/validation/build-teletubbies-Desktop_Qt_5_7_0_GCC_64bit-Release/data_1503097732603.txt"); n = 2^24 # N AR2 0.99
data = read.table("/home/marius/Dokumenter/master/blocking/validation/build-teletubbies-Desktop_Qt_5_7_0_GCC_64bit-Release/data_1503105060703.txt"); n = 2^20 # G AR1 0.99
data = read.table("/home/marius/Dokumenter/master/blocking/validation/build-teletubbies-Desktop_Qt_5_7_0_GCC_64bit-Release/data_1503105124934.txt"); n = 2^20 # G AR2 0.99
data = read.table("/home/marius/Dokumenter/master/blocking/validation/build-teletubbies-Desktop_Qt_5_7_0_GCC_64bit-Release/data_1503105875943.txt"); n = 2^27 # G AR1 0.99
data = read.table("/home/marius/Dokumenter/master/blocking/validation/build-teletubbies-Desktop_Qt_5_7_0_GCC_64bit-Release/AR2Log2n18_56472.txt"); n = 2^18 # G 0.99
data = read.table("/home/marius/Dokumenter/master/blocking/validation/build-teletubbies-Desktop_Qt_5_7_0_GCC_64bit-Release/AR1Log2n24_56245.txt"); n = 2^24 # G 0.99
data = read.table("/home/marius/Dokumenter/master/blocking/validation/build-teletubbies-Desktop_Qt_5_7_0_GCC_64bit-Release/AR2Log2n24_46310.txt"); n = 2^24 # N 0.99


new = (((data$V2-data$V1)/data$V1 )^2)
old = (((data$V3-data$V1)/data$V1 )^2)

#new[new>2] = NA
#old[old>2] = NA

phi = data$V4; explainatory = n/(1/-log(phi))  #ar1

#ar2
phi1 = data$V4; phi2 = data$V5; 
x = -phi1/(2*phi2); y = -(- phi1^2 - 4*phi2)^0.5/(2*phi2);
z = (-phi2)^(-0.5); theta = atan2(y,x); b = atan( (phi1/(z*(1 - phi2)) - cos(theta) )/sin(theta) )
a = 1/(cos(b))
explainatory = n*(-log(-phi2))/(2*(1 + log(abs(a) )))
#explainatory[explainatory < 2e2] = NA

# ---------


old.fit = glm(old ~ log(explainatory), family = Gamma(link = log) )
new.fit = glm(new ~ log(explainatory), family = Gamma(link = log) )

summary(old.fit)
summary(new.fit)

newdata = data.frame(explainatory = seq(min(explainatory, na.rm=T),max(explainatory, na.rm=T),max(explainatory, na.rm=T)/100000 ) )
new.pred = predict(new.fit, newdata=newdata,se.fit = T, type="response")
old.pred = predict(old.fit, newdata=newdata,se.fit = T, type="response")


plot(explainatory + runif(length(explainatory), 0,20),   old,pch=8,cex=0.01,col="red", xlab = "n/tau", log = "x", ylab = "rel.error square",xlim=c( 5e1,1e6) ,ylim=c(0.001, 1) )
points(explainatory + runif(length(explainatory),0,20), new,pch=8,cex=0.01,col="blue" )
lines(newdata$explainatory, new.pred$fit, cex=1, lwd=2, col="blue")
#lines(newdata$explainatory, (new.pred$fit + 1.96*new.pred$se.fit), lwd=2, lty=3, col="blue")
#lines(newdata$explainatory, (new.pred$fit - 1.96*new.pred$se.fit)^8, lwd=2, lty=3, col="blue")

lines(newdata$explainatory, old.pred$fit, cex=1, lwd=1.5, col="red")
#lines(newdata$explainatory, (pred2o$fit + 1.96*pred$se.fit), lwd=2, lty=3)
#lines(newdata$explainatory, (pred2o$fit - 1.96*pred$se.fit), lwd=2, lty=3)
