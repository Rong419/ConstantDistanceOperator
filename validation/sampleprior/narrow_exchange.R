args = commandArgs(trailingOnly=TRUE)

logfile.folder <- args[1]
output.figure.folder  <- args[2]

#logfile.folder <- "/Users/rzha419/Workspace/ConstantDistanceOperator/out/"
#output.figure.folder <- "/Users/rzha419/Workspace/ConstantDistanceOperator/validation/sampleprior/"

#This script is to validate the NarrowExchange operator on a tree with 3 taxa
#The distance among taxa are all equal to 10 
##############

#LogNormal distribution of rates
M <- 0.0
S <- 0.4

#the two constant rates
r1 <- 1.0
r2 <- 1.0

#the rate prior density to sample
densityFromRates <- function (r.p) {
  r.x <- (5 - r.p * 5) / 10;
  dlnorm(r1,meanlog=M, sdlog=S) * dlnorm(r2,meanlog=M, sdlog=S) * dlnorm(r.p,meanlog=M, sdlog=S) * dlnorm(r.x,meanlog=M, sdlog=S); 
}

#integrate the density function
#the area of the curve
A=integrate(densityFromRates,0,1)
B=A$value

#normalize the density function
#make sure that Normalized should be equal to 1
G <- function (r.p) {
  densityFromRates(r.p)/B; 
}
Normalized <- integrate(G,0,1)

#the mean of the distribution
E <- function (r.p) {
  r.p * G(r.p); 
}

mean = integrate(E,0,1)

#the standard deviation of the distribution
Mean = mean$value
F <- function (r.p) {
  (r.p - Mean) * (r.p - Mean) * G(r.p); 
}

variance = integrate(F, 0, 1)
Stdev = sqrt(variance$value)

l <- read.table(file=paste0(logfile.folder,"testNarrow.log"), sep="\t", header=T)

pdf (paste0(output.figure.folder,"NarrowExchange.pdf"),width=7,height=5)
hist(l$rates.4,prob=T,breaks=80,xlab=c("rp"),main='')
curve(G,from=0,to=1, col="red",add=TRUE,lwd=2)
dev.off()




