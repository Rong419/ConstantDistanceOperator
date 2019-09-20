library(ape)

args = commandArgs(trailingOnly=TRUE)

logfile.folder <- args[1]
output.figure.folder  <- args[2]

#logfile.folder <- "/Users/rzha419/Workspace/ConstantDistanceOperator/out/artifacts/ConstantDistanceOperator_jar/"
#output.figure.folder <- "/Users/rzha419/Workspace/ConstantDistanceOperator/out/artifacts/ConstantDistanceOperator_jar/"

#This script is to validate the operations of SimpleDistance
#that proposes a root time 
#and maintains genetic distances of the two branches linked to the root
##############

#initialize the constant distance and internal node time
d1 <- 0.1; d2 <- 0.2; d3 <- 0.4; d4 <- 0.27
t <- 1;

#functions of rates on the 4 branches
#i.e. distance divided by the time duration
r1 <- 0.1;
r2 <- 0.2;
r3 <- function (T) { d3 / T;}
r4 <- function (T) { d4 / (T-t);}

#LogNormal distribution of rates
M = -3;
S = 0.25;

#population size
N = 0.3

#distribution of rates
densityFromRates <- function (T) {
dlnorm(r1,meanlog=M, sdlog=S) * dlnorm(r2,meanlog=M, sdlog=S) * dlnorm(r3(T),meanlog=M, sdlog=S) * dlnorm(r4(T),meanlog=M, sdlog=S);
}

#distribution of times
densityFromTimes <- function (T) {
(1 / N) * exp(-3 * t / N) * (1 / N) * exp(-(T-t) / N);
}

#distribution of tree
densityFromTree <- function (T) {
densityFromTimes(T) * densityFromRates(T);
}

#integrate the density function
#the area of the curve
A=integrate(densityFromTree,1,15)
B=A$value

#normalize the density function
#make sure that Normalized should be equal to 1
G <- function (T) {
 	densityFromTree(T)/B; 
 }
Normalized <- integrate(G,1,15)

#the mean of the distribution
E <- function (T) {
	 T*G(T); 
 }
 
mean = integrate(E,1,15)

#the standard deviation of the distribution
Mean = mean$value
F <- function (T) {
  (T - Mean) * (T - Mean) * G(T); 
  }
  
variance = integrate(F, 1, 15)
Stdev = sqrt(variance$value)

Table = c(paste0("M = ", Mean," with absolute error < ", mean$abs.error, ", S = ", Stdev))

write.table(Table,file=paste0(output.figure.folder,"numerical_root_SD.txt"))


l <- read.table(file=paste0(logfile.folder,"root_SD.log"), sep="\t", header=T)

pdf (paste0(output.figure.folder,"root_SD.pdf"),width=7,height=5)
hist(l$TreeHeight,prob=T,breaks=80,xlab="",main='')
curve(G,from=1,to=15, col="red",add=TRUE,lwd=2)
dev.off()
