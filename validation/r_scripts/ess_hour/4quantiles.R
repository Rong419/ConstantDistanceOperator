#This R script is used to compare the ESS of the parameters
#Series 4
#Deal with the case where "quantiles" are used for rates

l <- read.table(file="/Users/rzha419/Desktop/posterior.txt", sep="\t", header=T);
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
posterior4<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/likelihood.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
likelihood4<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/prior.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
prior4<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/freq1.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
freq41<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/freq2.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
freq42<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/freq3.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
freq43<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/freq4.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
freq44<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/treeheight.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
treeheight4<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/ratemean.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
ratemean4<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/ratevariance.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
ratevariance4<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/ratecoef.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
ratecoef4<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/kapa.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
kapa4<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/birthrate.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
birthrate4<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/deathrate.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}  
deathrate4<-as.numeric(h)


M4=matrix(nrow = 14, ncol = 1)
M4[1]=mean(posterior4)
M4[2]=mean(likelihood4)
M4[3]=mean(prior4)
M4[4]=mean(treeheight4)
M4[5]=mean(ratemean4)
M4[6]=mean(ratevariance4)
M4[7]=mean(ratecoef4)
M4[8]=mean(kapa4)
M4[9]=mean(birthrate4)
M4[10]=mean(deathrate4)
M4[11]=mean(freq41)
M4[12]=mean(freq42)
M4[13]=mean(freq43)
M4[14]=mean(freq44)


