# This script is used to compare the estimated value(1) and the real value(2)
# (1) sampled parameters by using simulated sequences
# (2) samples obtained from prior distribution of the parameter



########################################
pdf ('120compare.pdf',width=10,height=6)
par(mfrow=c(2,3))
#get real values of the parameters that are used to simulate the sequences
TH <- read.table(file="~/Desktop/TreH.txt",sep="\t")
TH=as.matrix(TH)

#get estimated parameters from the MCMC samples by simulated sequences
H1 <- read.table(file="~/Desktop/TreeHeight.txt",sep="\t",header=T)

#get 95% HPD
K1=H1[8,]
K1=K1[2:101]
#get mean values
G1=H1[1,]
G1=G1[2:101]

#estimated mean
G=matrix(nrow=100,ncol=1)
for (i in 1:100) {
  x=factor(G1[1,i]);
  G[i]=levels(x);
}
G=as.numeric(G)

#plot real value vs estimated value
par(mar = c(5, 6, 3, 3))
plot(TH,G,xlim=c(1,5),ylim=c(1,5),pch=4,col='red',xlab='Real value',ylab='Estimated value',xaxt='n',yaxt='n',cex.lab=1.8)
axis(1,at=seq(1,5,1),labels=seq(1,5,1),cex.axis = 1.8)
axis(2,at=seq(1,5,1),labels=seq(1,5,1),cex.axis = 1.8)
title(main="TreeHeight",cex.main=1.8)

#plot 95%HPD
for (i in 1:100) {
  a=as.character(K1[1,i]);
  b=strsplit(a,",");
  c=substr(b[[1]][1],2,7);
  d=substr(b[[1]][2],2,7);
  c=as.numeric(c);
  d=as.numeric(d);
  x=c(TH[i],TH[i])
  y=c(c,d)
  par(new=TRUE)
  plot(x,y,pch=1,xlab='',ylab='',xlim=c(1,5),ylim=c(1,5),xaxt='n',yaxt='n')
  segments(TH[i],c,TH[i],d)
}
abline(0,1,col='blue',lwd=2)

#######################
TH <- read.table(file="~/Desktop/Kap.txt",sep="\t")
TH=as.matrix(TH)

H1 <- read.table(file="~/Desktop/Kappa.txt",sep="\t",header=T)
K1=H1[8,]
K1=K1[2:101]
G1=H1[1,]
G1=G1[2:101]

G=matrix(nrow=100,ncol=1)
for (i in 1:100) {
  x=factor(G1[1,i]);
  G[i]=levels(x);
}
G=as.numeric(G)
par(mar = c(5, 6, 3, 3))
plot(TH,G,xlim=c(0,1.5),ylim=c(0,1.5),pch=4,col='red',xlab='Real value',ylab='Estimated value',xaxt='n',yaxt='n',cex.lab=1.8)
axis(1,cex.axis = 1.8)
axis(2,cex.axis = 1.8)
title(main="Kappa",cex.main=1.8)

for (i in 1:100) {
  a=as.character(K1[1,i]);
  b=strsplit(a,",");
  c=substr(b[[1]][1],2,7);
  d=substr(b[[1]][2],2,7);
  c=as.numeric(c);
  d=as.numeric(d);
  x=c(TH[i],TH[i])
  y=c(c,d)
  par(new=TRUE)
  plot(x,y,pch=1,xlab='',ylab='',xlim=c(0,1.5),ylim=c(0,1.5),xaxt='n',yaxt='n')
  segments(TH[i],c,TH[i],d)
}
abline(0,1,col='blue',lwd=2)

###########################
TH <- read.table(file="~/Desktop/Bir.txt",sep="\t")
TH=as.matrix(TH)

H1 <- read.table(file="~/Desktop/BirthRate.txt",sep="\t",header=T)
K1=H1[8,]
K1=K1[2:101]
G1=H1[1,]
G1=G1[2:101]

G=matrix(nrow=100,ncol=1)
for (i in 1:100) {
  x=factor(G1[1,i]);
  G[i]=levels(x);
}
G=as.numeric(G)
par(mar = c(5, 6, 3, 3))
plot(TH,G,xlim=c(1,4),ylim=c(1,4),pch=4,col='red',xlab='Real value',ylab='Estimated value',xaxt='n',yaxt='n',cex.lab=1.8)
axis(1,at=seq(1,4,1),labels=seq(1,4,1),cex.axis = 1.8)
axis(2,at=seq(1,4,1),labels=seq(1,4,1),cex.axis = 1.8)
title(main="BirthRate",cex.main=1.8)

for (i in 1:100) {
  a=as.character(K1[1,i]);
  b=strsplit(a,",");
  c=substr(b[[1]][1],2,7);
  d=substr(b[[1]][2],2,7);
  c=as.numeric(c);
  d=as.numeric(d);
  x=c(TH[i],TH[i])
  y=c(c,d)
  par(new=TRUE)
  plot(x,y,pch=1,xlab='',ylab='',xlim=c(1,4),ylim=c(1,4),xaxt='n',yaxt='n')
  segments(TH[i],c,TH[i],d)
}
abline(0,1,col='blue',lwd=2)

######################
TH <- read.table(file="~/Desktop/RatM.txt",sep="\t")
TH=as.matrix(TH)

H1 <- read.table(file="~/Desktop/RateMean.txt",sep="\t",header=T)
K1=H1[8,]
K1=K1[2:101]
G1=H1[1,]
G1=G1[2:101]

G=matrix(nrow=100,ncol=1)
for (i in 1:100) {
  x=factor(G1[1,i]);
  G[i]=levels(x);
}
G=as.numeric(G)

par(mar = c(5, 6, 3, 3))
plot(0,0,xlim=c(0,100),ylim=c(0.9,1.2),xlab='Number of runs',ylab='Estimated value',cex.lab=1.8,yaxt='n',xaxt='n')
axis(2,at=seq(0.9,1.2,0.1),labels=seq(0.9,1.2,0.1),cex.axis = 1.8)
axis(1,at=seq(0,100,25),labels=seq(0,100,25),cex.axis = 1.8)
title(main="RateMean",cex.main=1.8)
#The following code is for experiments that don't sample rate mean
for (i in 1:100) {
  a=as.character(K1[1,i]);
  b=strsplit(a,",");
  c=substr(b[[1]][1],2,7);
  d=substr(b[[1]][2],2,7);
  c=as.numeric(c);
  d=as.numeric(d);
  x=c(i,i)
  y=c(c,d)
  par(new=TRUE)
  plot(x,y,pch=1,xlab='',ylab='',xlim=c(0,100),ylim=c(0.9,1.2),xaxt='n',yaxt='n')
  segments(i,c,i,d)
}
abline(h=1.0,col="red",lwd=2)


############################
TH <- read.table(file="~/Desktop/Ucld.txt",sep="\t")
TH=as.matrix(TH)

H1 <- read.table(file="~/Desktop/UcldStdev.txt",sep="\t",header=T)
K1=H1[8,]
K1=K1[2:101]
G1=H1[1,]
G1=G1[2:101]

G=matrix(nrow=100,ncol=1)
for (i in 1:100) {
  x=factor(G1[1,i]);
  G[i]=levels(x);
}
G=as.numeric(G)

par(mar = c(5, 6, 3, 3))
#plot(TH,G,xlim=c(0,0.6),ylim=c(0,0.6),pch=4,col='red',xlab='Real value',ylab='Eatimated value')
plot(TH,G,xlim=c(0,0.6),ylim=c(0,0.6),pch=4,col='red',xlab='Real value',ylab='Estimated value',xaxt='n',yaxt='n',cex.lab=1.8)
axis(1,at=seq(0,0.6,0.3),labels=seq(0,0.6,0.3),cex.axis = 1.8)
axis(2,at=seq(0,0.6,0.3),labels=seq(0,0.6,0.3),cex.axis = 1.8)
title(main="S3",cex.main=1.8)

for (i in 1:100) {
  a=as.character(K1[1,i]);
  b=strsplit(a,",");
  c=substr(b[[1]][1],2,7);
  d=substr(b[[1]][2],2,7);
  c=as.numeric(c);
  d=as.numeric(d);
  x=c(TH[i],TH[i])
  y=c(c,d)
  par(new=TRUE)
  plot(x,y,pch=1,xlab='',ylab='',xlim=c(0,0.6),ylim=c(0,0.6),xaxt='n',yaxt='n')
  segments(TH[i],c,TH[i],d)
}
abline(0,1,col='blue',lwd=2)

###########################
TH <- read.table(file="~/Desktop/Fre1.txt",sep="\t")
TH=as.matrix(TH)

H1 <- read.table(file="~/Desktop/Freq1.txt",sep="\t",header=T)
K1=H1[8,]
K1=K1[2:101]
G1=H1[1,]
G1=G1[2:101]

G=matrix(nrow=100,ncol=1)
for (i in 1:100) {
  x=factor(G1[1,i]);
  G[i]=levels(x);
}
G=as.numeric(G)

par(mar = c(5, 6, 3, 3))
plot(TH,G,xlim=c(0.1,0.5),ylim=c(0.1,0.5),pch=4,col='red',xlab='Real value',ylab='Estimated value',xaxt='n',yaxt='n',cex.lab=1.8)
axis(1,at=seq(0.1,0.5,0.2),labels=seq(0.1,0.5,0.2),cex.axis = 1.8)
axis(2,at=seq(0.1,0.5,0.2),labels=seq(0.1,0.5,0.2),cex.axis = 1.8)
title(main="Frequency",cex.main=1.8)

for (i in 1:100) {
  a=as.character(K1[1,i]);
  b=strsplit(a,",");
  c=substr(b[[1]][1],2,7);
  d=substr(b[[1]][2],2,7);
  c=as.numeric(c);
  d=as.numeric(d);
  x=c(TH[i],TH[i])
  y=c(c,d)
  par(new=TRUE)
  plot(x,y,pch=1,xlab='',ylab='',xlim=c(0.1,0.5),ylim=c(0.1,0.5),xaxt='n',yaxt='n')
  segments(TH[i],c,TH[i],d)
}
abline(0,1,col='blue',lwd=2)
dev.off()





#################################################################
pdf ('20compare.pdf',width=10,height=6)
par(mfrow=c(2,3))
#get real values of the parameters that are used to simulate the sequences
TH <- read.table(file="/Users/ryanzhang/Desktop/TreH.txt",sep="\t")
TH=as.matrix(TH)

#get estimated parameters from the MCMC samples by simulated sequences
H1 <- read.table(file="/Users/ryanzhang/Desktop/TreeHeight.txt",sep="\t",header=T)

#get 95% HPD
K1=H1[8,]
K1=K1[2:101]
#get mean values
G1=H1[1,]
G1=G1[2:101]

#estimated mean
G=matrix(nrow=100,ncol=1)
for (i in 1:100) {
  x=factor(G1[1,i]);
  G[i]=levels(x);
}
G=as.numeric(G)

#plot real value vs estimated value
par(mar = c(5, 6, 3, 3))
plot(TH,G,xlim=c(0,3),ylim=c(0,3),pch=4,col='red',xlab='Real value',ylab='Estimated value',xaxt='n',yaxt='n',cex.lab=1.8)
axis(1,at=seq(0,3,1),labels=seq(0,3,1),cex.axis = 1.8)
axis(2,at=seq(0,3,1),labels=seq(0,3,1),cex.axis = 1.8)
title(main="TreeHeight",cex.main = 1.8)

#plot 95%HPD
for (i in 1:100) {
  a=as.character(K1[1,i]);
  b=strsplit(a,",");
  c=substr(b[[1]][1],2,7);
  d=substr(b[[1]][2],2,7);
  c=as.numeric(c);
  d=as.numeric(d);
  x=c(TH[i],TH[i])
  y=c(c,d)
  par(new=TRUE)
  plot(x,y,pch=1,xlab='',ylab='',xlim=c(0,3),ylim=c(0,3),xaxt='n',yaxt='n')
  segments(TH[i],c,TH[i],d)
}
abline(0,1,col='blue',lwd=2)

#######################
TH <- read.table(file="/Users/ryanzhang/Desktop/Kap.txt",sep="\t")
TH=as.matrix(TH)

H1 <- read.table(file="/Users/ryanzhang/Desktop/Kappa.txt",sep="\t",header=T)
K1=H1[8,]
K1=K1[2:101]
G1=H1[1,]
G1=G1[2:101]

G=matrix(nrow=100,ncol=1)
for (i in 1:100) {
  x=factor(G1[1,i]);
  G[i]=levels(x);
}
G=as.numeric(G)
par(mar = c(5, 6, 3, 3))
plot(TH,G,xlim=c(0,1.5),ylim=c(0,1.5),pch=4,col='red',xlab='Real value',ylab='Estimated value',xaxt='n',yaxt='n',cex.lab=1.8)
axis(1,cex.axis = 1.8)
axis(2,cex.axis = 1.8)
title(main="Kappa",cex.main = 1.8)

for (i in 1:100) {
  a=as.character(K1[1,i]);
  b=strsplit(a,",");
  c=substr(b[[1]][1],2,7);
  d=substr(b[[1]][2],2,7);
  c=as.numeric(c);
  d=as.numeric(d);
  x=c(TH[i],TH[i])
  y=c(c,d)
  par(new=TRUE)
  plot(x,y,pch=1,xlab='',ylab='',xlim=c(0,1.5),ylim=c(0,1.5),xaxt='n',yaxt='n')
  segments(TH[i],c,TH[i],d)
}
abline(0,1,col='blue',lwd=2)

###########################
TH <- read.table(file="/Users/ryanzhang/Desktop/Bir.txt",sep="\t")
TH=as.matrix(TH)

H1 <- read.table(file="/Users/ryanzhang/Desktop/BirthRate.txt",sep="\t",header=T)
K1=H1[8,]
K1=K1[2:101]
G1=H1[1,]
G1=G1[2:101]

G=matrix(nrow=100,ncol=1)
for (i in 1:100) {
  x=factor(G1[1,i]);
  G[i]=levels(x);
}
G=as.numeric(G)
par(mar = c(5, 6, 3, 3))
plot(TH,G,xlim=c(1,4),ylim=c(1,4),pch=4,col='red',xlab='Real value',ylab='Estimated value',xaxt='n',yaxt='n',cex.lab=1.8)
axis(1,at=seq(1,4,1),labels=seq(1,4,1),cex.axis = 1.8)
axis(2,at=seq(1,4,1),labels=seq(1,4,1),cex.axis = 1.8)
title(main="BirthRate",cex.main = 1.8)

for (i in 1:100) {
  a=as.character(K1[1,i]);
  b=strsplit(a,",");
  c=substr(b[[1]][1],2,7);
  d=substr(b[[1]][2],2,7);
  c=as.numeric(c);
  d=as.numeric(d);
  x=c(TH[i],TH[i])
  y=c(c,d)
  par(new=TRUE)
  plot(x,y,pch=1,xlab='',ylab='',xlim=c(1,4),ylim=c(1,4),xaxt='n',yaxt='n')
  segments(TH[i],c,TH[i],d)
}
abline(0,1,col='blue',lwd=2)

######################
TH <- read.table(file="/Users/ryanzhang/Desktop/RatM.txt",sep="\t")
TH=as.matrix(TH)

H1 <- read.table(file="/Users/ryanzhang/Desktop/RateMean.txt",sep="\t",header=T)
K1=H1[8,]
K1=K1[2:101]
G1=H1[1,]
G1=G1[2:101]

G=matrix(nrow=100,ncol=1)
for (i in 1:100) {
  x=factor(G1[1,i]);
  G[i]=levels(x);
}
G=as.numeric(G)
par(mar = c(5, 6, 3, 3))
plot(0,0,xlim=c(0,100),ylim=c(0.85,1.25),xlab='Number of runs',ylab='Estimated value',cex.lab=1.8,yaxt='n',xaxt='n')
axis(2,at=seq(0.85,1.25,0.15),labels=seq(0.85,1.25,0.15),cex.axis = 1.8)
axis(1,at=seq(0,100,25),labels=seq(0,100,25),cex.axis = 1.8)
title(main="RateMean",cex.main = 1.8)
#The following code is for experiments that don't sample rate mean
for (i in 1:100) {
  a=as.character(K1[1,i]);
  b=strsplit(a,",");
  c=substr(b[[1]][1],2,7);
  d=substr(b[[1]][2],2,7);
  c=as.numeric(c);
  d=as.numeric(d);
  x=c(i,i)
  y=c(c,d)
  par(new=TRUE)
  plot(x,y,pch=1,xlab='',ylab='',xlim=c(0,100),ylim=c(0.85,1.25),xaxt='n',yaxt='n')
  segments(i,c,i,d)
}
#title(main="RateMean",xlab='Number of runs',ylab='Estimated value')
abline(h=1.0,col="red",lwd=2)


############################
TH <- read.table(file="/Users/ryanzhang/Desktop/Ucld.txt",sep="\t")
TH=as.matrix(TH)

H1 <- read.table(file="/Users/ryanzhang/Desktop/UcldStdev.txt",sep="\t",header=T)
K1=H1[8,]
K1=K1[2:101]
G1=H1[1,]
G1=G1[2:101]

G=matrix(nrow=100,ncol=1)
for (i in 1:100) {
  x=factor(G1[1,i]);
  G[i]=levels(x);
}
G=as.numeric(G)
par(mar = c(5, 6, 3, 3))
#plot(TH,G,xlim=c(0,0.6),ylim=c(0,0.6),pch=4,col='red',xlab='Real value',ylab='Eatimated value')
plot(TH,G,xlim=c(0,0.6),ylim=c(0,0.6),pch=4,col='red',xlab='Real value',ylab='Estimated value',xaxt='n',yaxt='n',cex.lab=1.8)
axis(1,at=seq(0,0.6,0.3),labels=seq(0,0.6,0.3),cex.axis = 1.8)
axis(2,at=seq(0,0.6,0.3),labels=seq(0,0.6,0.3),cex.axis = 1.8)
title(main="S3",cex.main = 1.8)

for (i in 1:100) {
  a=as.character(K1[1,i]);
  b=strsplit(a,",");
  c=substr(b[[1]][1],2,7);
  d=substr(b[[1]][2],2,7);
  c=as.numeric(c);
  d=as.numeric(d);
  x=c(TH[i],TH[i])
  y=c(c,d)
  par(new=TRUE)
  plot(x,y,pch=1,xlab='',ylab='',xlim=c(0,0.6),ylim=c(0,0.6),xaxt='n',yaxt='n')
  segments(TH[i],c,TH[i],d)
}
abline(0,1,col='blue',lwd=2)

###########################
TH <- read.table(file="/Users/ryanzhang/Desktop/Fre1.txt",sep="\t")
TH=as.matrix(TH)

H1 <- read.table(file="/Users/ryanzhang/Desktop/Freq1.txt",sep="\t",header=T)
K1=H1[8,]
K1=K1[2:101]
G1=H1[1,]
G1=G1[2:101]

G=matrix(nrow=100,ncol=1)
for (i in 1:100) {
  x=factor(G1[1,i]);
  G[i]=levels(x);
}
G=as.numeric(G)
par(mar = c(5, 6, 3, 3))
plot(TH,G,xlim=c(0.1,0.5),ylim=c(0.1,0.5),pch=4,col='red',xlab='Real value',ylab='Estimated value',xaxt='n',yaxt='n',cex.lab=1.8)
axis(1,at=seq(0.1,0.5,0.2),labels=seq(0.1,0.5,0.2),cex.axis = 1.8)
axis(2,at=seq(0.1,0.5,0.2),labels=seq(0.1,0.5,0.2),cex.axis = 1.8)
title(main="Frequency",cex.main = 1.8)

for (i in 1:100) {
  a=as.character(K1[1,i]);
  b=strsplit(a,",");
  c=substr(b[[1]][1],2,7);
  d=substr(b[[1]][2],2,7);
  c=as.numeric(c);
  d=as.numeric(d);
  x=c(TH[i],TH[i])
  y=c(c,d)
  par(new=TRUE)
  plot(x,y,pch=1,xlab='',ylab='',xlim=c(0.1,0.5),ylim=c(0.1,0.5),xaxt='n',yaxt='n')
  segments(TH[i],c,TH[i],d)
}
abline(0,1,col='blue',lwd=2)
dev.off()




#################################################################
TH <- read.table(file="/Users/ryanzhang/Desktop/TreL.txt",sep="\t")
TH=as.matrix(TH)

H1 <- read.table(file="/Users/ryanzhang/Desktop/TreeLength.txt",sep="\t",header=T)
K1=H1[8,]
K1=K1[2:101]
G1=H1[1,]
G1=G1[2:101]

G=matrix(nrow=100,ncol=1)
for (i in 1:100) {
  x=factor(G1[1,i]);
  G[i]=levels(x);
}
G=as.numeric(G)

plot(TH,G,xlim=c(0,20),ylim=c(0,20),pch=4,col='red',xlab='Real value',ylab='Eatimated value')
title(main="TreeLength")
for (i in 1:100) {
  a=as.character(K1[1,i]);
  b=strsplit(a,",");
  c=substr(b[[1]][1],2,7);
  d=substr(b[[1]][2],2,7);
  c=as.numeric(c);
  d=as.numeric(d);
  x=c(TH[i],TH[i])
  y=c(c,d)
  par(new=TRUE)
  plot(x,y,pch=1,xlab='',ylab='',xlim=c(0,20),ylim=c(0,20))
  segments(TH[i],c,TH[i],d)
}
abline(0,1,col='blue',lwd=2)


plot(TH,G,xlim=c(0.9,1.1),ylim=c(0.9,1.1),pch=4,col='red',xlab='Real value',ylab='Eatimated value')
title(main="RateMean",cex.main = 1.6)

for (i in 1:100) {
  a=as.character(K1[1,i]);
  b=strsplit(a,",");
  c=substr(b[[1]][1],2,7);
  d=substr(b[[1]][2],2,7);
  c=as.numeric(c);
  d=as.numeric(d);
  x=c(TH[i],TH[i])
  y=c(c,d)
  par(new=TRUE)
  plot(x,y,pch=1,xlab='',ylab='',xlim=c(0.9,1.1),ylim=c(0.9,1.1))
  segments(TH[i],c,TH[i],d)
}
abline(0,1,col='blue',lwd=2)
#the line of y=0+1*x (real value equals to estimated value)
abline(0,1,col='blue',lwd=2)