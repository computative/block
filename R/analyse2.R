data = read.table("/home/marius/Dokumenter/master/blocking/validation/build-teletubbies-Desktop_Qt_5_7_0_GCC_64bit-Release/data_1501957507286.txt")
data = read.table("/home/marius/Dokumenter/master/blocking/validation/build-teletubbies-Desktop_Qt_5_7_0_GCC_64bit-Release/data_1502047901965.txt")
data = read.table("/home/marius/Dokumenter/master/blocking/validation/build-teletubbies-Desktop_Qt_5_7_0_GCC_64bit-Release/data_1502834807876.txt")
data = read.table("/home/marius/Dokumenter/master/blocking/validation/build-teletubbies-Desktop_Qt_5_7_0_GCC_64bit-Release/data_1502834846949.txt")
data = read.table("/home/marius/Dokumenter/master/blocking/validation/build-teletubbies-Desktop_Qt_5_7_0_GCC_64bit-Release/data_1502885662776.txt")
new = (((data$V2-data$V1)/data$V1)^2)
old = (((data$V3-data$V1)/data$V1)^2)

new[new > 10] = NA
old[old > 10] = NA

explainatory = 2^26/(1/-log(data$V4))  #ar1
explainatory = 2^26/(2/-log(-data$V5)) #ar2

# head

new = new[order(explainatory)]
old = old[order(explainatory)]
explainatory = explainatory[order(explainatory)]
explainatory = explainatory[1:440]
new = new[1:440]
old = old[1:440]

#tail

new = new[order(explainatory)]
old = old[order(explainatory)]
explainatory = explainatory[order(explainatory)]
explainatory = explainatory[440:3808]
new = new[440:3808]
old = old[440:3808]




old.fit = glm(old ~ log(explainatory), family = poisson(link=log))
new.fit = glm(new ~ log(explainatory), family = poisson(link=log))
summary(old.fit)
summary(new.fit)

newdata = data.frame(explainatory = seq(0,max(explainatory),1000) )
new.pred = predict(new.fit, newdata=newdata,se.fit = T, type="response")
old.pred = predict(old.fit, newdata=newdata,se.fit = T, type="response")

plot(explainatory, new,pch=8,cex=0.05,col="blue", xlab = "obs/tau", ylab = "rel.error square", log="x", ylim=c(0,0.6) )
points(explainatory, old,pch=8,cex=0.01,col="red" )
lines(newdata$explainatory, new.pred$fit,cex=1, lwd=1.5, col="blue")
#lines(newdata$explainatory, new.pred$fit + 1.96*new.pred$se.fit, lwd=2, lty=3)
#lines(newdata$explainatory, new.pred$fit - 1.96*new.pred$se.fit, lwd=2, lty=3)

lines(newdata$explainatory, old.pred$fit,cex=1, lwd=1.5, col="red")
#lines(newdata$explainatory, (old.pred$fit + 1.96*old.pred$se.fit), lwd=2, lty=3)
#lines(newdata$explainatory, (old.pred$fit - 1.96*old.pred$se.fit), lwd=2, lty=3)