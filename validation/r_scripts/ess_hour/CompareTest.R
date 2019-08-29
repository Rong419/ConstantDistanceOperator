pdf ('Efficiency.pdf',width=16,height=10)
par(mfrow=c(2,2))
###########################################################################
##20 taxa
T1 <- read.table(file="~/Desktop/Error_categories1.txt",sep="\t")
T2 <- read.table(file="~/Desktop/Error_categories2.txt",sep="\t")
T3 <- read.table(file="~/Desktop/Error_categories3.txt",sep="\t")
T1=as.matrix(T1[1])
T2=as.matrix(T2[1])
T3=as.matrix(T3[1])

S1 <- read.table(file="~/Desktop/Error_cons1.txt",sep="\t")
S2 <- read.table(file="~/Desktop/Error_cons2.txt",sep="\t")
S3 <- read.table(file="~/Desktop/Error_cons3.txt",sep="\t")
S1=as.matrix(S1[1])
S2=as.matrix(S2[1])
S3=as.matrix(S3[1])

P1 <- read.table(file="~/Desktop/Error_nocons1.txt",sep="\t")
P2 <- read.table(file="~/Desktop/Error_nocons2.txt",sep="\t")
P3 <- read.table(file="~/Desktop/Error_nocons3.txt",sep="\t")
P1=as.matrix(P1[1])
P2=as.matrix(P2[1])
P3=as.matrix(P3[1])

par(mar = c(3, 6, 3, 3))
boxplot(T1,ylim=c(4000,25000),xlim=c(1,3),boxwex=0.3,at=1.1,outline=FALSE,yaxt='n')
boxplot(T2,ylim=c(4000,25000),col="green",add=TRUE,xlim=c(1,3),boxwex=0.3,at=1.3,outline=FALSE,yaxt='n')
boxplot(T3,ylim=c(4000,25000),col="blue",add=TRUE,xlim=c(1,3),boxwex=0.3,at=1.5,outline=FALSE,yaxt='n')

boxplot(S1,ylim=c(4000,25000),add=TRUE,xlim=c(1,3),boxwex=0.3,at=1.8,outline=FALSE,yaxt='n')
boxplot(S2,ylim=c(4000,25000),col="green",add=TRUE,xlim=c(1,3),boxwex=0.3,at=2.0,outline=FALSE,yaxt='n')
boxplot(S3,ylim=c(4000,25000),col="blue",add=TRUE,xlim=c(1,3),boxwex=0.3,at=2.2,outline=FALSE,yaxt='n')

boxplot(P1,ylim=c(4000,25000),xlim=c(1,3),add=TRUE,boxwex=0.3,at=2.5,outline=FALSE,yaxt='n')
boxplot(P2,ylim=c(4000,25000),col="green",add=TRUE,xlim=c(1,3),boxwex=0.3,at=2.7,outline=FALSE,yaxt='n')
boxplot(P3,ylim=c(4000,25000),col="blue",add=TRUE,xlim=c(1,3),boxwex=0.3,at=2.9,outline=FALSE,yaxt='n')

abline(v=1.65,col="orange")
abline(v=2.35,col="orange")

axis(1,at=c(1.3,2,2.7),labels=c("categories","cons","nocons"),cex.axis = 2.2)
title(main="Runing time of MCMC (taxa=20)",ylab ="time (second)",cex.main=2.2,cex.lab=2.2)
axis(2,cex.axis = 2.2)
#legend("topright",inset=.01,c("long","medium","short"),pch=1,col=c(1,3,4),box.lty = 0,cex=1.5,bg='transparent')
#B=c(mean(T1),mean(T2),mean(T3),mean(S1),mean(S2),mean(S3),mean(P1),mean(P2),mean(P3))
#write.csv(B,file="/Users/rzha419/Desktop/b.csv")


L1 <- read.table(file="~/Desktop/UcldStdev_categories1.txt", sep="\t", header=T)
L2 <- read.table(file="~/Desktop/UcldStdev_categories2.txt", sep="\t", header=T)
L3 <- read.table(file="/Users/ryanzhang/Desktop/UcldStdev_categories3.txt", sep="\t", header=T)

Q1 <- read.table(file="~/Desktop/UcldStdev_cons1.txt", sep="\t", header=T)
Q2 <- read.table(file="~/Desktop/UcldStdev_cons2.txt", sep="\t", header=T)
Q3 <- read.table(file="~/Desktop/UcldStdev_cons3.txt", sep="\t", header=T)

R1 <- read.table(file="~/Desktop/UcldStdev_nocons1.txt", sep="\t", header=T)
R2 <- read.table(file="~/Desktop/UcldStdev_nocons2.txt", sep="\t", header=T)
R3 <- read.table(file="/Users/ryanzhang/Desktop/UcldStdev_nocons3.txt", sep="\t", header=T)

L1=L1[10,]
L2=L2[10,]
L3=L3[10,]

Q1=Q1[10,]
Q2=Q2[10,]
Q3=Q3[10,]

R1=R1[10,]
R2=R2[10,]
R3=R3[10,]

h1 = matrix(nrow = 50, ncol = 1);
for (i in 1:50){
  X=L1[1,i+1];
  Y=factor(X);
  h1[i]=levels(Y);
}
UcldStdev1<-as.numeric(h1)

h2 = matrix(nrow = 50, ncol = 1);
for (i in 1:50){
  X=L2[1,i+1];
  Y=factor(X);
  h2[i]=levels(Y);
}
UcldStdev2<-as.numeric(h2)

h3 = matrix(nrow = 50, ncol = 1);
for (i in 1:50){
  X=L3[1,i+1];
  Y=factor(X);
  h3[i]=levels(Y);
}
UcldStdev3<-as.numeric(h3)

h4 = matrix(nrow = 50, ncol = 1);
for (i in 1:50){
  X=Q1[1,i+1];
  Y=factor(X);
  h4[i]=levels(Y);
}
UcldStdev4<-as.numeric(h4)

h5 = matrix(nrow = 50, ncol = 1);
for (i in 1:50){
  X=Q2[1,i+1];
  Y=factor(X);
  h5[i]=levels(Y);
}
UcldStdev5<-as.numeric(h5)

h6 = matrix(nrow = 50, ncol = 1);
for (i in 1:50){
  X=Q3[1,i+1];
  Y=factor(X);
  h6[i]=levels(Y);
}
UcldStdev6<-as.numeric(h6)

h7 = matrix(nrow = 50, ncol = 1);
for (i in 1:50){
  X=R1[1,i+1];
  Y=factor(X);
  h7[i]=levels(Y);
}
UcldStdev7<-as.numeric(h7)

h8 = matrix(nrow = 50, ncol = 1);
for (i in 1:50){
  X=R2[1,i+1];
  Y=factor(X);
  h8[i]=levels(Y);
}
UcldStdev8<-as.numeric(h8)

h9 = matrix(nrow = 50, ncol = 1);
for (i in 1:50){
  X=R3[1,i+1];
  Y=factor(X);
  h9[i]=levels(Y);
}
UcldStdev9<-as.numeric(h9)
par(mar = c(3, 6, 3, 3))
boxplot(UcldStdev1,ylim=c(1,1300),xlim=c(1,3),boxwex=0.3,at=1.1,outline=FALSE,yaxt='n')
boxplot(UcldStdev2,ylim=c(1,1300),col="green",add=TRUE,xlim=c(1,3),boxwex=0.3,at=1.3,outline=FALSE,yaxt='n')
boxplot(UcldStdev3,ylim=c(1,1300),col="blue",add=TRUE,xlim=c(1,3),boxwex=0.3,at=1.5,outline=FALSE,yaxt='n')

boxplot(UcldStdev4,ylim=c(1,1300),xlim=c(1,3),add=TRUE,boxwex=0.3,at=1.8,outline=FALSE,yaxt='n')
boxplot(UcldStdev5,ylim=c(1,1300),col="green",add=TRUE,xlim=c(1,3),boxwex=0.3,at=2,outline=FALSE,yaxt='n')
boxplot(UcldStdev6,ylim=c(1,1300),col="blue",add=TRUE,xlim=c(1,3),boxwex=0.3,at=2.2,outline=FALSE,yaxt='n')

boxplot(UcldStdev7,ylim=c(1,1300),xlim=c(1,3),add=TRUE,boxwex=0.3,at=2.5,outline=FALSE,yaxt='n')
boxplot(UcldStdev8,ylim=c(1,1300),col="green",add=TRUE,xlim=c(1,3),boxwex=0.3,at=2.7,outline=FALSE,yaxt='n')
boxplot(UcldStdev9,ylim=c(1,1300),col="blue",add=TRUE,xlim=c(1,3),boxwex=0.3,at=2.9,outline=FALSE,yaxt='n')

abline(v=1.65,col="orange")
abline(v=2.35,col="orange")

axis(1,at=c(1.3,2,2.7),labels=c("categories","cons","nocons"),cex.axis = 2.2)
title(main="Ess of S3 (taxa=20)",ylab="ESS",cex.main=2.2,cex.lab=2.2)
axis(2,cex.axis = 2.2)



###########################################################################
##120 taxa
T1 <- read.table(file="~/Desktop/Error_categories4.txt",sep="\t")
T2 <- read.table(file="~/Desktop/Error_categories5.txt",sep="\t")
T3 <- read.table(file="~/Desktop/Error_categories6.txt",sep="\t")
T1=as.matrix(T1[1])
T2=as.matrix(T2[1])
T3=as.matrix(T3[1])

S1 <- read.table(file="~/Desktop/Error_cons4.txt",sep="\t")
S2 <- read.table(file="~/Desktop/Error_cons5.txt",sep="\t")
S3 <- read.table(file="~/Desktop/Error_cons6.txt",sep="\t")
S1=as.matrix(S1[1])
S2=as.matrix(S2[1])
S3=as.matrix(S3[1])

P1 <- read.table(file="~/Desktop/Error_nocons4.txt",sep="\t")
P2 <- read.table(file="~/Desktop/Error_nocons5.txt",sep="\t")
P3 <- read.table(file="~/Desktop/Error_nocons6.txt",sep="\t")
P1=as.matrix(P1[1])
P2=as.matrix(P2[1])
P3=as.matrix(P3[1])

par(mar = c(3, 6, 3, 3))
boxplot(T1,ylim=c(10000,60000),xlim=c(1,3),boxwex=0.3,at=1.1,outline=FALSE,yaxt='n')
boxplot(T2,ylim=c(10000,60000),col="green",add=TRUE,xlim=c(1,3),boxwex=0.3,at=1.3,outline=FALSE,yaxt='n')
boxplot(T3,ylim=c(10000,60000),col="blue",add=TRUE,xlim=c(1,3),boxwex=0.3,at=1.5,outline=FALSE,yaxt='n')

boxplot(S1,ylim=c(10000,60000),add=TRUE,xlim=c(1,3),boxwex=0.3,at=1.8,outline=FALSE,yaxt='n')
boxplot(S2,ylim=c(10000,60000),col="green",add=TRUE,xlim=c(1,3),boxwex=0.3,at=2.0,outline=FALSE,yaxt='n')
boxplot(S3,ylim=c(10000,60000),col="blue",add=TRUE,xlim=c(1,3),boxwex=0.3,at=2.2,outline=FALSE,yaxt='n')

boxplot(P1,ylim=c(10000,60000),xlim=c(1,3),add=TRUE,boxwex=0.3,at=2.5,outline=FALSE,yaxt='n')
boxplot(P2,ylim=c(10000,60000),col="green",add=TRUE,xlim=c(1,3),boxwex=0.3,at=2.7,outline=FALSE,yaxt='n')
boxplot(P3,ylim=c(10000,60000),col="blue",add=TRUE,xlim=c(1,3),boxwex=0.3,at=2.9,outline=FALSE,yaxt='n')

abline(v=1.65,col="orange")
abline(v=2.35,col="orange")

axis(1,at=c(1.3,2,2.7),labels=c("categories","cons","nocons"),cex.axis = 2.2)
title(main="Runing time of MCMC (taxa=120)",ylab ="time (second)",cex.main=2.2,cex.lab=2.2)
axis(2,cex.axis = 2.2)
#legend("topright",inset=.01,c("long","medium","short"),pch=1,col=c(1,3,4),box.lty = 0,cex=1.5,bg='transparent')
#B=c(mean(T1),mean(T2),mean(T3),mean(S1),mean(S2),mean(S3),mean(P1),mean(P2),mean(P3))
#write.csv(B,file="/Users/rzha419/Desktop/b.csv")


L1 <- read.table(file="~/Desktop/UcldStdev_categories4.txt", sep="\t", header=T)
L2 <- read.table(file="~/Desktop/UcldStdev_categories5.txt", sep="\t", header=T)
L3 <- read.table(file="~/Desktop/UcldStdev_categories6.txt", sep="\t", header=T)

Q1 <- read.table(file="~/Desktop/UcldStdev_cons4.txt", sep="\t", header=T)
Q2 <- read.table(file="~/Desktop/UcldStdev_cons5.txt", sep="\t", header=T)
Q3 <- read.table(file="~/Desktop/UcldStdev_cons6.txt", sep="\t", header=T)

R1 <- read.table(file="~/Desktop/UcldStdev_nocons4.txt", sep="\t", header=T)
R2 <- read.table(file="~/Desktop/UcldStdev_nocons5.txt", sep="\t", header=T)
R3 <- read.table(file="~/Desktop/UcldStdev_nocons6.txt", sep="\t", header=T)

L1=L1[10,]
L2=L2[10,]
L3=L3[10,]

Q1=Q1[10,]
Q2=Q2[10,]
Q3=Q3[10,]

R1=R1[10,]
R2=R2[10,]
R3=R3[10,]

h1 = matrix(nrow = 50, ncol = 1);
for (i in 1:50){
  X=L1[1,i+1];
  Y=factor(X);
  h1[i]=levels(Y);
}
UcldStdev1<-as.numeric(h1)

h2 = matrix(nrow = 50, ncol = 1);
for (i in 1:50){
  X=L2[1,i+1];
  Y=factor(X);
  h2[i]=levels(Y);
}
UcldStdev2<-as.numeric(h2)

h3 = matrix(nrow = 50, ncol = 1);
for (i in 1:50){
  X=L3[1,i+1];
  Y=factor(X);
  h3[i]=levels(Y);
}
UcldStdev3<-as.numeric(h3)

h4 = matrix(nrow = 50, ncol = 1);
for (i in 1:50){
  X=Q1[1,i+1];
  Y=factor(X);
  h4[i]=levels(Y);
}
UcldStdev4<-as.numeric(h4)

h5 = matrix(nrow = 50, ncol = 1);
for (i in 1:50){
  X=Q2[1,i+1];
  Y=factor(X);
  h5[i]=levels(Y);
}
UcldStdev5<-as.numeric(h5)

h6 = matrix(nrow = 50, ncol = 1);
for (i in 1:50){
  X=Q3[1,i+1];
  Y=factor(X);
  h6[i]=levels(Y);
}
UcldStdev6<-as.numeric(h6)

h7 = matrix(nrow = 50, ncol = 1);
for (i in 1:50){
  X=R1[1,i+1];
  Y=factor(X);
  h7[i]=levels(Y);
}
UcldStdev7<-as.numeric(h7)

h8 = matrix(nrow = 50, ncol = 1);
for (i in 1:50){
  X=R2[1,i+1];
  Y=factor(X);
  h8[i]=levels(Y);
}
UcldStdev8<-as.numeric(h8)

h9 = matrix(nrow = 50, ncol = 1);
for (i in 1:50){
  X=R3[1,i+1];
  Y=factor(X);
  h9[i]=levels(Y);
}
UcldStdev9<-as.numeric(h9)
par(mar = c(3, 6, 3, 3))
boxplot(UcldStdev1,ylim=c(0,350),xlim=c(1,3),boxwex=0.3,at=1.1,outline=FALSE,yaxt='n')
boxplot(UcldStdev2,ylim=c(0,350),col="green",add=TRUE,xlim=c(1,3),boxwex=0.3,at=1.3,outline=FALSE,yaxt='n')
boxplot(UcldStdev3,ylim=c(0,350),col="blue",add=TRUE,xlim=c(1,3),boxwex=0.3,at=1.5,outline=FALSE,yaxt='n')

boxplot(UcldStdev4,ylim=c(0,350),xlim=c(1,3),add=TRUE,boxwex=0.3,at=1.8,outline=FALSE,yaxt='n')
boxplot(UcldStdev5,ylim=c(0,350),col="green",add=TRUE,xlim=c(1,3),boxwex=0.3,at=2,outline=FALSE,yaxt='n')
boxplot(UcldStdev6,ylim=c(0,350),col="blue",add=TRUE,xlim=c(1,3),boxwex=0.3,at=2.2,outline=FALSE,yaxt='n')

boxplot(UcldStdev7,ylim=c(0,350),xlim=c(1,3),add=TRUE,boxwex=0.3,at=2.5,outline=FALSE,yaxt='n')
boxplot(UcldStdev8,ylim=c(0,350),col="green",add=TRUE,xlim=c(1,3),boxwex=0.3,at=2.7,outline=FALSE,yaxt='n')
boxplot(UcldStdev9,ylim=c(0,350),col="blue",add=TRUE,xlim=c(1,3),boxwex=0.3,at=2.9,outline=FALSE,yaxt='n')

abline(v=1.65,col="orange")
abline(v=2.35,col="orange")

axis(1,at=c(1.3,2,2.7),labels=c("categories","cons","nocons"),cex.axis = 2.2)
title(main="Ess of S3 (taxa=120)",ylab="ESS",cex.main=2.2,cex.lab=2.2)
axis(2,cex.axis = 2.2)
#legend("topright",inset=.01,c("long","medium","short"),pch=1,col=c(1,3,4),box.lty = 0,cex=1.5,bg='transparent')

#M=c(mean(UcldStdev1),mean(UcldStdev2),mean(UcldStdev3),mean(UcldStdev4),mean(UcldStdev5),mean(UcldStdev6),mean(UcldStdev7),mean(UcldStdev8),mean(UcldStdev9))
#write.csv(M,file="/Users/rzha419/Desktop/a.csv")
dev.off()

################################################################
###########################################################################
##Real data set (primates)
T1 <- read.table(file="/Users/ryanzhang/Desktop/Error_categories.txt",sep="\t")

T2 <- read.table(file="/Users/ryanzhang/Desktop/Error_cons.txt",sep="\t")

T3 <- read.table(file="/Users/ryanzhang/Desktop/Error_nocons.txt",sep="\t")

T1=as.matrix(T1[1])
T2=as.matrix(T2[1])
T3=as.matrix(T3[1])

par(mar = c(5, 6, 3, 3))
boxplot(T1,ylim=c(10000,25000),xlim=c(1,3),boxwex=0.3,at=1+0.5,outline=FALSE,xaxt="n",yaxt='n')
boxplot(T2,ylim=c(10000,25000),col="green",add=TRUE,xlim=c(1,3),boxwex=0.3,at=2,outline=FALSE,xaxt="n",yaxt='n')
boxplot(T3,ylim=c(10000,25000),col="blue",add=TRUE,xlim=c(1,3),boxwex=0.3,at=3-0.5,outline=FALSE,xaxt="n",yaxt='n')
axis(1,at=c(1.5,2,2.5),labels=c("categories","cons","nocons"),cex.axis=1.8)
title(main="Runing time of MCMC (83 primates)",ylab="time (second)",cex.lab=1.8,cex.main=1.8)
axis(2,at=seq(10000,25000,5000),labels=seq(10000,25000,5000),cex.axis = 1.8)
A=c(1.5,2,2.5)
#B=c(mean(T1),mean(T2),mean(T3))
B=c(18368,15046,18404)
par(new=TRUE)
plot(A,B,ylim=c(10000,25000),xlim=c(1,3),type='o',xlab='',ylab='',col="orange",xaxt="n",yaxt='n')

text(A,B,labels=B,pos=3,cex=1.8)



L1 <- read.table(file="/Users/ryanzhang/Desktop/UcldStdev_categories.txt", sep="\t", header=T)
L2 <- read.table(file="/Users/ryanzhang/Desktop/UcldStdev_cons.txt", sep="\t", header=T)
L3 <- read.table(file="/Users/ryanzhang/Desktop/UcldStdev_nocons.txt", sep="\t", header=T)

C=L1[10,]
D=L2[10,]
E=L3[10,]

h1 = matrix(nrow = 30, ncol = 1);
for (i in 1:50){
  X=C[1,i+1];
  Y=factor(X);
  h1[i]=levels(Y);
}
UcldStdev1<-as.numeric(h1)

h2 = matrix(nrow = 30, ncol = 1);
for (i in 1:50){
  X=D[1,i+1];
  Y=factor(X);
  h2[i]=levels(Y);
}
UcldStdev2<-as.numeric(h2)

h3 = matrix(nrow = 30, ncol = 1);
for (i in 1:50){
  X=E[1,i+1];
  Y=factor(X);
  h3[i]=levels(Y);
}
UcldStdev3<-as.numeric(h3)
par(mar = c(5, 6, 3, 3))
boxplot(UcldStdev1,ylim=c(0,300),xlim=c(1,3),boxwex=0.3,at=1+0.5,outline=FALSE,xaxt="n",yaxt='n')
boxplot(UcldStdev2,ylim=c(0,300),col="green",add=TRUE,xlim=c(1,3),boxwex=0.3,at=2,outline=FALSE,xaxt="n",yaxt='n')
boxplot(UcldStdev3,ylim=c(0,300),col="blue",add=TRUE,xlim=c(1,3),boxwex=0.3,at=3-0.5,outline=FALSE,xaxt="n",yaxt='n')
axis(1,at=c(1.5,2,2.5),labels=c("categories","cons","nocons"),cex.axis=1.8)
axis(2,at=seq(0,300,100),labels=seq(0,300,100),cex.axis = 1.8)
title(main="Ess of S3 (83 primates)",ylab="ESS",cex.lab=1.8,cex.main=1.8)

A=c(1.5,2,2.5)
#M=c(mean(UcldStdev1),mean(UcldStdev2),mean(UcldStdev3))

par(new=TRUE)
plot(A,M,ylim=c(0,300),xlim=c(1,3),type='o',xlab='',ylab='',col="orange",xaxt="n",yaxt='n')
M=c(9,133,15)
text(A,M,labels=M,pos=3,cex=1.8)



