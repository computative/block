mydata = read.table("/home/marius/Dokumenter/fys4411/postmann.pat/code/dataget/data/data.txt") 
a = mydata$V1

mydata = read.table("/home/marius/Dokumenter/fys4411/postmann.pat/code/dataget/data/data.txt") 
b = mydata$V1

mydata = read.table("/home/marius/Dokumenter/fys4411/postmann.pat/code/dataget/data/data.txt") 
c = mydata$V1

mydata = read.table("/home/marius/Dokumenter/fys4411/postmann.pat/code/dataget/data/data.txt") 
d = mydata$V1

mydata = read.table("/home/marius/Dokumenter/fys4411/postmann.pat/code/dataget/data/data.txt") 
e = mydata$V1

mydata = read.table("/home/marius/Dokumenter/fys4411/postmann.pat/code/dataget/data/data.txt") 
f = mydata$V1


mydata = read.table("/home/marius/Dokumenter/fys4411/postmann.pat/code/dataget/data/data.txt") 
g = mydata$V1

# yet todo
mydata = read.table("/home/marius/Dokumenter/fys4411/postmann.pat/code/dataget/data/data.txt") 
h = mydata$V1

library(boot)

X[2] = var(tsboot(g,mean,R=2^4, sim="geom", l=2^11)$t)
acf(g,15000)