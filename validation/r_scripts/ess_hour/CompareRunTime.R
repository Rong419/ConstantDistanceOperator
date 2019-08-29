T1 <- read.table(file="/Users/rzha419/Desktop/Out_categories.txt",sep="\t")

T2 <- read.table(file="/Users/rzha419/Desktop/Out_cons.txt",sep="\t")

T3 <- read.table(file="/Users/rzha419/Desktop/Out_nocons.txt",sep="\t")

T4 <- read.table(file="/Users/rzha419/Desktop/Out_quantiles.txt",sep="\t")

T1=as.matrix(T1[1])
T2=as.matrix(T2[1])
T3=as.matrix(T3[1])
T3=as.matrix(T4[1])

boxplot(T1,ylim=c(5000,8000),xlim=c(1,3),boxwex=0.3,at=1+0.5,outline=FALSE)
boxplot(T2,ylim=c(5000,8000),col="green",add=TRUE,xlim=c(1,3),boxwex=0.3,at=2,outline=FALSE)
boxplot(T3,ylim=c(5000,8000),col="blue",add=TRUE,xlim=c(1,3),boxwex=0.3,at=3-0.5,outline=FALSE)
axis(1,at=c(1.5,2,2.5),labels=c("categories","cons","nocons"))
title(main="Runing time of a short chain (length=5000)",ylab="time (second)")

A=c(1.5,2,2.5)
B=c(mean(T1),mean(T2),mean(T3))
par(new=TRUE)
plot(A,B,ylim=c(5000,8000),xlim=c(1,3),type='o',xlab='',ylab='',col="orange",xaxt="n",yaxt='n')



plot(T1,type='o',xlab = "run",ylim=c(5400,8000), ylab = "time",xaxt="n")
par(new=TRUE)
plot(T2,type='o',xlab = "",ylim=c(5400,8000),col="red", ylab = "",xaxt="n")
par(new=TRUE)
plot(T3,type='o',xlab = "",ylim=c(5400,8000),col="green", ylab = "",xaxt="n")