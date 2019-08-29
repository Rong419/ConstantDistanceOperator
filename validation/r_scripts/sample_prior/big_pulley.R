#This R script is used to validate the operations of BigPulley
#Test the tree with 3 taxa, i.e. A B and C
#There are 3 topologies in total
#The integration includes 3 parts for corresponding 3 topologies

###############################################
##1. only consider the tree ((C B) A)

#Initialize distance
d1 <- 1; d2 <- 2; 
#d3 <- 0.4; d4 <- 0.2;

#all the possible values for genetic distance d4
Y = seq(0,d1,length.out=1000);

##time of the higher child of the root
t <- 5;

##the constant popolation size in Coalescent model
N = 0.3;

##the lognormal distribution that rates follow
M = -3;
S = 0.25;

##root time T
X = seq(0,5,length.out=1000);

##initialize the parameters
P=matrix(nrow=1000,ncol=1);
P3=matrix(nrow=1000,ncol=1);
P4=matrix(nrow=1000,ncol=1);


##integrate d4 out
j = 1;
for (T in X) {
			i = 1;
			r2 = d2 / T;
			r1 = (0.4 + 0.2) / T;
			t1 = (1 / N) * exp(-3 * T / N);
		    t2 = (1 / N) * exp(-(t-T) / N);
		    P2 = log(t1) + log(t2);
		for (d4 in Y){
			r3 = (d1 - d4) / t;
		    #r2 = (d2 - d4) / t;
		    r4 = d4 / (t - T);
		    P1 = log(dlnorm(r1,meanlog=M, sdlog=S)) + log(dlnorm(r2,meanlog=M, sdlog=S)) + log(dlnorm(r3,meanlog=M, sdlog=S)) + log(dlnorm(r4,meanlog=M, sdlog=S)); 
		    LogP = P1 + P2;
		    P[i] = exp(LogP);
		    i = i + 1;
		}
		P3[j]=mean(P);
        P4[j]=P3[j]*T;
        j=j+1;	  
}
Num=seq(1,999,length.out=999)
A=matrix(nrow=999,ncol=1);
B=matrix(nrow=999,ncol=1);
C=matrix(nrow=999,ncol=1);
k=1;
for (x in Num) {
	A[k] = P3[x]; #Mean(LogP)
	B[k] = X[x];#T
	C[k] = P4[x];#T*Mean(LogP)
	k = k + 1;
}
plot(B,A);

#Normalize constant
Cons = sum(A);
#mean
E = sum(C/Cons);
#deviation
e = (B - E) * (B - E);
V = sum(e * A / Cons);
#standard deviation
std = sqrt(V);
##Validation:D=1
D = sum(A / Cons);

#integrate T out
for (d4 in Y) {
	        i = 1;
	        r3 = (d1 - d4) / t;   
   for (T in X) { 
		    r2 = d2 / T;
		    r1 = (0.4 + 0.2) / T; 
		    r4 = d4 / (t - T);	 
		    t1 = (1 / N) * exp(-3 * T / N);
		    t2 = (1 / N) * exp(-(t-T) / N);
		    P1 = log(dlnorm(r1,meanlog=M, sdlog=S)) + log(dlnorm(r2,meanlog=M, sdlog=S)) + log(dlnorm(r3,meanlog=M, sdlog=S)) + log(dlnorm(r4,meanlog=M, sdlog=S));
		    P2 = log(t1) + log(t2);
		    LogP = P1 + P2;
		    P[i] = exp(LogP);
		    i = i + 1;
		}
		P3[j]=mean(P);
        P4[j]=P3[j]*d4;
        j=j+1;	  
}
Num=seq(2,1000,length.out=999)
A=matrix(nrow=999,ncol=1);
B=matrix(nrow=999,ncol=1);
C=matrix(nrow=999,ncol=1);
k=1;
for (x in Num) {
	A[k] = P3[x]; #Mean(LogP)
	B[k] = Y[x];#d4
	C[k] = P4[x];#d4*Mean(LogP)
	k = k + 1;
}
plot(B,A);

#Normalize constant
Cons = sum(A);
#mean
E = sum(C/Cons);
#variance
e = (B - E) * (B - E);
V = sum(e * A / Cons);
#standard devation
std = sqrt(V);
##Validation:D=1
D = sum(A / Cons);

###############################################
#2. only consider the big pulley (Tree: (A C) B)
##Initialize distance
d1 <- 1; d2 <- 2; 
#d3 <- 0.4; d4 <- 0.2;

##genetic distance d4
Y = seq(0,d2,length.out=1000);

##time of the higher child of the root
t <- 5;

##the constant popolation size Coalescene model
N = 0.3;

##the lognormal distribution that rates follow
M = -3;
S = 0.25;

##root time T
X = seq(0,5,length.out=1000);

##initialize the parameters
P=matrix(nrow=1000,ncol=1);
P3=matrix(nrow=1000,ncol=1);
P4=matrix(nrow=1000,ncol=1);
j = 1;
##integrate d4
for (T in X) {
			i = 1;
			r1 = d1 / T;
			r2 = (0.4 + 0.2) / T;
			t1 = (1 / N) * exp(-3 * T / N);
		    t2 = (1 / N) * exp(-(t-T) / N);
		    P2 = log(t1) + log(t2);
		for (d4 in Y){
		    r3 = (d2 - d4) / t;
		    r4 = d4 / (t - T);
		    P1 = log(dlnorm(r1,meanlog=M, sdlog=S)) + log(dlnorm(r2,meanlog=M, sdlog=S)) + log(dlnorm(r3,meanlog=M, sdlog=S)) + log(dlnorm(r4,meanlog=M, sdlog=S)); 
		    LogP = P1 + P2;
		    P[i] = exp(LogP);
		    i = i + 1;
		}
		P3[j]=mean(P);
        P4[j]=P3[j]*T;
        j=j+1;	  
}
Num=seq(1,999,length.out=999)
A=matrix(nrow=999,ncol=1);
B=matrix(nrow=999,ncol=1);
C=matrix(nrow=999,ncol=1);
k=1;
for (x in Num) {
	A[k] = P3[x]; #Mean(LogP)
	B[k] = X[x];#T
	C[k] = P4[x];#T*Mean(LogP)
	k = k + 1;
}
plot(B,A);

#Normalize constant
Cons = sum(A);
#mean
E = sum(C/Cons);
#variance
e = (B - E) * (B - E);
V = sum(e * A / Cons);
#standard devation
std = sqrt(V);
##Validation:D=1
D = sum(A / Cons);


##read files from MCMC, .log and .trees
library(ape)

strees <- read.nexus("/Users/rzha419/Documents/SoftWare/beast2-master/dna1.subst.trees")

e <- sapply(strees, function(x) {x$edge.length[1]})

hist(e,prob=T, breaks=100)

l <- read.table(file="/Users/rzha419/Documents/SoftWare/beast2-master/dna.log", sep="\t", header=T)

hist(l$mrcatime.ingroup.,prob=T,breaks=80)

