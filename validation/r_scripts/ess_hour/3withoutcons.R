#This R script is used to compare the ESS of the parameters
#Series 3
#Deal with the case where rates are used directly
#NO constant distance operators are included

l <- read.table(file="/Users/rzha419/Desktop/posterior.txt", sep="\t", header=T);
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
posterior3<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/likelihood.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
likelihood3<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/prior.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
prior3<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/freq1.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
freq31<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/freq2.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
freq32<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/freq3.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
freq33<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/freq4.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
freq34<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/treeheight.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
treeheight3<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/ratemean.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
ratemean3<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/ratevariance.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
ratevariance3<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/ratecoef.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
ratecoef3<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/kappa.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
kapa3<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/birthrate.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
birthrate3<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/deathrate.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}  
deathrate3<-as.numeric(h)


M3=matrix(nrow = 14, ncol = 1)
M3[1]=mean(posterior3)
M3[2]=mean(likelihood3)
M3[3]=mean(prior3)
M3[4]=mean(treeheight3)
M3[5]=mean(ratemean3)
M3[6]=mean(ratevariance3)
M3[7]=mean(ratecoef3)
M3[8]=mean(kapa3)
M3[9]=mean(birthrate3)
M3[10]=mean(deathrate3)
M3[11]=mean(freq31)
M3[12]=mean(freq32)
M3[13]=mean(freq33)
M3[14]=mean(freq34)




