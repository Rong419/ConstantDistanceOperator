#This R script is used to compare the ESS of the parameters
#Series 1
#Deal with the case where rates are used directly
#And the constant distance operators are included


l <- read.table(file="/Users/rzha419/Desktop/posterior.txt", sep="\t", header=T);
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
posterior1<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/likelihood.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
likelihood1<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/prior.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
prior1<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/freq1.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
freq11<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/freq2.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
freq12<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/freq3.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
freq13<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/freq4.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
freq14<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/treeheight.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
treeheight1<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/ratemean.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
ratemean1<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/ratevariance.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
ratevariance1<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/ratecoef.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
ratecoef1<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/kappa.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
kapa1<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/birthrate.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
birthrate1<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/deathrate.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}  
deathrate1<-as.numeric(h)

M1=matrix(nrow = 14, ncol = 1);
M1[1]=mean(posterior1)
M1[2]=mean(likelihood1)
M1[3]=mean(prior1)
M1[4]=mean(treeheight1)
M1[5]=mean(ratemean1)
M1[6]=mean(ratevariance1)
M1[7]=mean(ratecoef1)
M1[8]=mean(kapa1)
M1[9]=mean(birthrate1)
M1[10]=mean(deathrate1)
M1[11]=mean(freq11)
M1[12]=mean(freq12)
M1[13]=mean(freq13)
M1[14]=mean(freq14)


