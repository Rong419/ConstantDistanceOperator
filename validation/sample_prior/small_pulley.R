library(ape)

args = commandArgs(trailingOnly=TRUE)

substitution.tree.folder <- args[1]
output.figure.folder <- args[2]

#substitution.tree.folder <- "/Users/rzha419/Workspace/ConstantDistanceOperator/out/artifacts/ConstantDistanceOperator_jar/"
#output.figure.folder <- "/Users/rzha419/Workspace/ConstantDistanceOperator/out/artifacts/ConstantDistanceOperator_jar/"

#This R script is used to validate the operations of SmallPulley
#that proposes one genetic distance 
#and maintains genetic distances of the two branches linked to the root
####################################################
#integrate d3 out
#initalization
#the distance d1,d2 and the sum of d3 and d4 (D=d3+d4)
d1 <- 0.1; d2 <- 0.2; D <- 0.67;

#root time
T <- 10;

#the node time of the child of the root
t <- 1;

#the rates on branches
r1 <- 0.1;
r2 <- 0.2;
r3 <- function (d3) { d3 / T;}
r4 <- function (d3) { (D-d3) / (T-t);}

#Lognormal distribution for rates
M = -3;
S = 0.25;

#population size
N = 0.3

#distribution of rates
densityFromRates <- function (d3) {
dlnorm(r1,meanlog=M, sdlog=S) * dlnorm(r2,meanlog=M, sdlog=S) * dlnorm(r3(d3),meanlog=M, sdlog=S) * dlnorm(r4(d3),meanlog=M, sdlog=S);
}

#distribution of times
densityFromTimes = (1 / N) * exp(-3 * t / N) * (1 / N) * exp(-(T-t) / N);


#distribution of tree
densityFromTree <- function (d3) {
densityFromTimes * densityFromRates(d3);
}

#integrate the density function
#the area of the curve
A=integrate(densityFromTree,0,D)
B=A$value

#normalize the density function
#make sure that Normalized should be equal to 1
G<- function (d3) {
 	densityFromTree(d3)/B; 
 }
 
Normalized = integrate(G,0,D)

#the mean of the distribution
E<- function (d3) {
	 d3*G(d3); 
 }
 
mean = integrate(E,0,D)

#the standard deviation of the distribution
Mean = mean$value
F <- function (d3) {
  (d3-Mean)*(d3-Mean)*G(d3); 
 }

variance = integrate(F,0,D)
Stdev = sqrt(variance$value)

Table = c(paste0("SmallPulley: M = ", Mean, "with absolute error < ", mean$abs.error, ", S = ", Stdev))

write.table(Table,file=paste0(output.figure.folder,"numerical_root_SP.txt"))


####################################################
#read tree files
strees <- read.nexus(paste0(substitution.tree.folder,"subst_root_SP.trees"))
#get the distance on branch
e <- sapply(strees, function(x) {x$edge.length[4]})


pdf (paste0(output.figure.folder,"root_SP.pdf"),width=7,height=5)
#plot the distribution of the distance
hist(e,prob=T, breaks=100,xlab="",main='')
#plot the integrated distribution
curve(G,from=0,to=D, col="red",add=TRUE)
dev.off()
