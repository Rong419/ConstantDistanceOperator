#This R script is used to compare the ESS of the parameters
#Series 2
#Deal with the case where "categoires" are used for rates


#read the log file
l <- read.table(file="/Users/rzha419/Desktop/posterior.txt", sep="\t", header=T);
#get the ESS line 
A=l[10,];
#save the values of ESS to "posterior2" in numeric format
#matrix "h" contains ESS in char format
#"nrow" equals to the number of ESS
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
posterior2<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/likelihood.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
likelihood2<-as.numeric(h)

################
l <- read.table(file="/Users/rzha419/Desktop/prior.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
prior2<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/freq1.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
freq21<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/freq2.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
freq22<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/freq3.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
freq23<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/freq4.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
freq24<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/treeheight.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
treeheight2<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/ratemean.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
ratemean2<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/ratevariance.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
ratevariance2<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/ratecoef.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
ratecoef2<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/kappa.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
kapa2<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/birthrate.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
birthrate2<-as.numeric(h)

l <- read.table(file="/Users/rzha419/Desktop/deathrate.txt", sep="\t", header=T)
A=l[10,];
h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}  
deathrate2<-as.numeric(h)

#get the mean of ESS for each parameter
#"nrow" of M2 is the number of parameters
M2=matrix(nrow = 14, ncol = 1)
M2[1]=mean(posterior2)
M2[2]=mean(likelihood2)
M2[3]=mean(prior2)
M2[4]=mean(treeheight2)
M2[5]=mean(ratemean2)
M2[6]=mean(ratevariance2)
M2[7]=mean(ratecoef2)
M2[8]=mean(kapa2)
M2[9]=mean(birthrate2)
M2[10]=mean(deathrate2)
M2[11]=mean(freq21)
M2[12]=mean(freq22)
M2[13]=mean(freq23)
M2[14]=mean(freq24)










