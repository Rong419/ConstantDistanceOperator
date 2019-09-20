library(ape)

logfile.folder <- args[1]
output.figure.folder  <- args[2]

#logfile.folder <- "/Users/rzha419/Workspace/ConstantDistanceOperator/out/artifacts/ConstantDistanceOperator_jar/" 
#output.figure.folder <- "/Users/rzha419/Workspace/ConstantDistanceOperator/out/artifacts/ConstantDistanceOperator_jar/" 


#This R script is used to validate the operations of BigPulley
#Test the tree with 3 taxa, i.e. A B and C
#There are 3 topologies in total
#The integration includes 3 parts for corresponding 3 topologies

###############################################
# Initilaize the rates and times
# r1 = 0.1; r2 = 0.1; r3 = 0.03; r4 = 0.04
# t = 5; T = 10;
# Calculate the distances
# d1 = r1 * t;
# d2 = r2 * t;
# d3 = r3 * T;
# d4 = r4 * (T - t);

d1 = 0.5; d2= 0.5; d3 = 0.3; d4 = 0.2;

N = 0.3; #population size
M = -3; #mean
S = 0.25; #standard deviation

Y = seq(0.001,0.499,length.out=50); # d4 range
Z = seq(0.1,10, length.out=50) # T range
X = seq(0,9.999,length.out=10000); # t range
P = matrix(nrow = 50, ncol = 1);
P1 = matrix(nrow = 10000, ncol = 1);
P2 = matrix(nrow = 10000, ncol = 1);
P3 = matrix(nrow = 50, ncol = 1);
P4 = matrix(nrow = 50, ncol = 1);
P5 = matrix(nrow = 50, ncol = 1);
P6 = matrix(nrow = 10000, ncol = 1);
P7 = matrix(nrow = 10000, ncol = 1);
P8 = matrix(nrow = 50, ncol = 1);
P9 = matrix(nrow = 50, ncol = 1);
P10 = matrix(nrow = 50, ncol = 1);
P11 = matrix(nrow = 50, ncol = 1);
P12 = matrix(nrow = 50, ncol = 1);
P13 = matrix(nrow = 50, ncol = 1);

k = 1;	
# iterating root time T
for (T_ in Z) { 
     j = 1;
     # iterating internal node time t
 for (t_ in X) {
	    i = 1;
	    a = (1 / N) * exp(-3 * t_ / N); # P(t1)
	    b = (1 / N) * exp(-(T_ - t_ ) / N); # P(t2)
	        # iterating distance d4
	        for (d4_ in Y) {
		       r1_ = d1 / t_;
		       r2_ = (d2 - d4_) / T_;
		       r3_ = (d3 + d4) / t_;
		       r4_ = d4_ / (T_ - t_);
		       P[i] =  a * b * dlnorm(r1_,meanlog=M, sdlog=S) * dlnorm(r2_,meanlog=M, sdlog=S) * dlnorm(r3_,meanlog=M, sdlog=S) * dlnorm(r4_,meanlog=M, sdlog=S);
		       P5[i] = P[i] * d4_;
               i = i + 1;
	           }	
	      P1[j] = mean(P5);
          P2[j] = P1[j] * t_;
          P6[j] = mean(P);
          P7[j] = P6[j] * t_;
          j = j + 1;	      
     }	
    P3[k] = mean(P1); ##Mean(Density)
    P4[k] = P3[k] * T_; ##T*Mean(Density)
    P10[k] = mean(P6);
    P11[k] = P10[k] * T_;

    P8[k] = mean(P7);
    P9[k] = P8[k] * T_;
    P12[k] = mean(P2);
    P13[k] = P12[k] * T_;
    k = k + 1;
}


plot(Z,P12);

##Normalize constant
C1=sum(P12); 
## mean of distribution
E1=sum(P13/C1);
##Variance and standard deviation
##B is the random variable X
e1=(Z-E1)*(Z-E1);
V1=sum(e1*P12/C1);
Sd1=sqrt(V1);
##Validation:D=1
D1=sum(P3/C1);

##(2) Integrate T and t out
d1 = 0.5; d2= 0.5; d3 = 0.3; d4 = 0.2;

N = 0.3; #population size
M = -3; #mean
S = 0.25; #standard deviation

Y = seq(0.001,0.499,length.out=50); # d4 range
Z = seq(1,10, length.out=50) # T range
P = matrix(nrow = 10000, ncol = 1);
P1 = matrix(nrow = 50, ncol = 1);
P2 = matrix(nrow = 50, ncol = 1);
P3 = matrix(nrow = 50, ncol = 1);
P4 = matrix(nrow = 50, ncol = 1);
P5 = matrix(nrow = 10000, ncol = 1);
P6 = matrix(nrow = 50, ncol = 1);
P7 = matrix(nrow = 50, ncol = 1);
P8 = matrix(nrow = 50, ncol = 1);
P9 = matrix(nrow = 50, ncol = 1);
P10 = matrix(nrow = 50, ncol = 1);
P11 = matrix(nrow = 50, ncol = 1);
P12 = matrix(nrow = 50, ncol = 1);
P13 = matrix(nrow = 50, ncol = 1);

k = 1;	
for (d4_ in Y) { 
     j = 1;
     for (T_ in Z) {
          X = seq(0,T_-0.01,length.out=10000); # t range
	    i = 1; 
	     for (t_ in X) {
                   a = (1 / N) * exp(-3 * t_ / N); # P(t1)
	             b = (1 / N) * exp(-(T_ - t_ ) / N); # P(t2)
		       r1_ = d1 / t_;
		       r2_ = (d2 - d4_) / T_;
		       r3_ = (d3 + d4) / t_;
		       r4_ = d4_ / (T_ - t_);
		       P[i] =  a * b * dlnorm(r1_,meanlog=M, sdlog=S) * dlnorm(r2_,meanlog=M, sdlog=S) * dlnorm(r3_,meanlog=M, sdlog=S) * dlnorm(r4_,meanlog=M, sdlog=S);
		       P5[i] = P[i] * t_;
                   i = i + 1;
	           }	
	    P1[j] = mean(P5);
          P2[j] = P1[j] * T_;
          P6[j] = mean(P);
          P7[j] = P6[j] * T_;
          j = j + 1;	      
     }	
    P3[k] = mean(P1); ##Mean(Density)
    P4[k] = P3[k] * d4_; ##T*Mean(Density)
    P10[k] = mean(P6);
    P11[k] = P10[k] * d4_;
    P8[k] = mean(P7);
    P9[k] = P8[k] * d4_;
    P12[k] = mean(P2);
    P13[k] = P12[k] * d4_;

    k = k + 1;
}

plot(Y,P12);

##Normalize constant
C1=sum(P12); 
## mean of distribution
E1=sum(P13/C1);
##Variance and standard deviation
##B is the random variable X
e1=(Y-E1)*(Y-E1);
V1=sum(e1*P10/C1);
Sd1=sqrt(V1);
##Validation:D=1
D1=sum(P3/C1);




##read files from MCMC, .log and .trees
strees <- read.nexus(paste0(logfile.folder,"subst_root_BP.trees"))
e <- sapply(strees, function(x) {x$edge.length[1]})
hist(e,prob=T, breaks=100,xlab="",main='')


l <- read.table(file=paste0(logfile.folder,"root_BP.log"), sep="\t", header=T)
hist(l$TreeHeight,prob=T,breaks=80)

