# read data
l <- read.table(file="/Users/rzha419/Desktop/Untitled.txt", sep="\t", header=T)

#get ESS for all the parameters
A=l[10,]

#delete the n/a ESS
B=A[,-c(8)]

#get the first parameter
x1=B[1,2]

#get ESS of the parameter
e1=levels(x1)[3]
#get the mean of the parameter
m1=levels(x1)[5]
#get the stderr of the mean of the parameter
s1=levels(x1)[2]
#get the std deviation of the parameter
d1=levels(x1)[4]
#B=A[,-c(8)];
B=A[,-c(8)];
#Step1: get the ESS
l6 <- read.table(file="/Users/rzha419/Desktop/6.txt", sep="\t", header=T);
A=l6[10,];

h = matrix(nrow = 10, ncol = 1);
for (i in 1:10){
  X=A[1,i+1];
  Y=factor(X);
  h[i]=levels(Y);
}
par(new=TRUE)
plot(e,xlab = "",type='o', ylab = "",ylim=c(0,6000),col="red")
plot(e,type='o',xlab = "parameter", ylab = "ESS",ylim=c(0,5000))


#Step2: get the rate mean 
M6 <- read.table(file="/Users/rzha419/Desktop/dna.log", sep="\t", header=T)
R6 = M6$rate.mean

#Step3: get the tree height 
H6 = M6$TreeHeight.dna
