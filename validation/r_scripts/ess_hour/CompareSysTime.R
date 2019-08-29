T1 <- read.table(file="/Users/rzha419/Desktop/Error_categories.txt",sep="\t")

T2 <- read.table(file="/Users/rzha419/Desktop/Error_cons.txt",sep="\t")

T3 <- read.table(file="/Users/rzha419/Desktop/Error_nocons.txt",sep="\t")

T1=as.matrix(T1[1])
T2=as.matrix(T2[1])
T3=as.matrix(T3[1])

boxplot(T1,ylim=c(15000,25000),xlim=c(1,3),boxwex=0.3,at=1+0.5,outline=FALSE)
boxplot(T2,ylim=c(15000,25000),col="green",add=TRUE,xlim=c(1,3),boxwex=0.3,at=2,outline=FALSE)
boxplot(T3,ylim=c(15000,25000),col="blue",add=TRUE,xlim=c(1,3),boxwex=0.3,at=3-0.5,outline=FALSE)
axis(1,at=c(1.5,2,2.5),labels=c("categories","cons","nocons"),cex.lab=1.8)
title(main="Runing time of MCMC (length=20000)",ylab="time (second)")

A=c(1.5,2,2.5)
B=c(mean(T1),mean(T2),mean(T3))
par(new=TRUE)
plot(A,B,ylim=c(15000,25000),xlim=c(1,3),type='o',xlab='',ylab='',col="orange",xaxt="n",yaxt='n')

text(A,B,labels=B,pos=2,cex=1)

L1 <- read.table(file="/Users/rzha419/Desktop/UcldStdev_categories.txt", sep="\t", header=T)
L2 <- read.table(file="/Users/rzha419/Desktop/UcldStdev_cons.txt", sep="\t", header=T)
L3 <- read.table(file="/Users/rzha419/Desktop/UcldStdev_nocons.txt", sep="\t", header=T)

C=L1[10,]
D=L2[10,]
E=L3[10,]

h1 = matrix(nrow = 50, ncol = 1);
for (i in 1:50){
  X=C[1,i+1];
  Y=factor(X);
  h1[i]=levels(Y);
}
UcldStdev1<-as.numeric(h1)

h2 = matrix(nrow = 50, ncol = 1);
for (i in 1:50){
  X=D[1,i+1];
  Y=factor(X);
  h2[i]=levels(Y);
}
UcldStdev2<-as.numeric(h2)

h3 = matrix(nrow = 50, ncol = 1);
for (i in 1:50){
  X=E[1,i+1];
  Y=factor(X);
  h3[i]=levels(Y);
}
UcldStdev3<-as.numeric(h3)

boxplot(UcldStdev1,ylim=c(1,1000),xlim=c(1,3),boxwex=0.3,at=1+0.5,outline=FALSE)
boxplot(UcldStdev2,ylim=c(1,1000),col="green",add=TRUE,xlim=c(1,3),boxwex=0.3,at=2,outline=FALSE)
boxplot(UcldStdev3,ylim=c(1,1000),col="blue",add=TRUE,xlim=c(1,3),boxwex=0.3,at=3-0.5,outline=FALSE)
axis(1,at=c(1.5,2,2.5),labels=c("categories","cons","nocons"))
title(main="Ess of ucldStdev (length=20000)",ylab="ESS")

A=c(1.5,2,2.5)
M=c(mean(UcldStdev1),mean(UcldStdev2),mean(UcldStdev3))

par(new=TRUE)
plot(A,M,ylim=c(1,1000),xlim=c(1,3),type='o',xlab='',ylab='',col="orange",xaxt="n",yaxt='n')

text(A,M,labels=M,pos=2,cex=1)
