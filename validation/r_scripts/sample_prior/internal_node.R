#This script is to validate the operations on internal nodes by
#1) reading the log file from MCMC, 
#2) presenting the histogram of the parameter
##############

#The tree has 3 taxa
#4 rates,4 distance
#2 node times (internal node: t,the root: T)


#initialize the constant distance and root time
d1 <- 0.1; d2 <- 0.2; d3 <- 0.4; d4 <- 0.27
T <- 10;

#functions of rates on the 4 branches
#i.e. distance divided by the time duration
r1 <- function (t) { d1 / t;}
r2 <- function (t) {d2 / t;}
r3 <- 0.04
r4 <- function (t) { d4 / (T-t);}

#LogNormal distribution of rates
M = -3;
S = 0.25;

#distribution of rates
#CDF in log space 
logRateDensity <- function (t) {
logd <- log(dlnorm(r1(t),meanlog=M, sdlog=S)) +log(dlnorm(r2(t),meanlog=M, sdlog=S)) +log(dlnorm(r3,meanlog=M, sdlog=S))+log(dlnorm(r4(t),meanlog=M, sdlog=S)) ;
}

#CDF in real space
densityFromRates <- function (t) {
exp(logRateDensity(t)); 
}

#integrate the density function
#the area of the curve
A=integrate(densityFromRates,0,T)
B=A$value

#normalize the density function
#each CDF divided by the whole area
G<- function (t) {
 	densityFromRates(t)/B; 
 }
integrate(G,0,T)

#the mean of the distribution
#i.e. m
E<- function (t) {
	 t*G(t); 
 }
 
integrate(E,0,T)

#the standard deviation of the distribution
# i.e. stdev
t1<- function (t) {
  (t-3.465785)*(t-3.465785); 
 }
 
F<- function (t) {
  t1(t)*G(t); 
  }

 integrate(F,0,T)
 
#read log file
l <- read.table(file="/Users/rzha419/desktop/dna.log", sep="\t", header=T)
#plot the histogram of the parameter
#i.e. the distribution of sampled root time 
 hist(l$mrcatime.ingroup.,prob=T,breaks=80,xlab="",main='')

#add the curve of the theoretical distribution of the root time
#the normalized CDF distribution
curve(G,from=0,to=T, col="red",add=TRUE,lwd=2)
curve(densityFromRates,from=0,to=T, col="red",add=TRUE)

d <- density(l$mrcatime.ingroup., n=length(l$mrcatime.ingroup.))
plot(d, main="")


################### Example ###################
d1 <- 0.4; d2 <- 0.8; d3 <- 2.4; d4 <- 1.6
 T <- 0.8;
 
 r1 <- function (t) { d1 / t;}
 r2 <- function (t) {d2 / t;}
 r3 <- 0.04
 r4 <- function (t) { d4 / (T-t);}
 
 M = -3;
 S = 0.25;
 
 logRateDensity <- function (t) {
   logd <- log(dlnorm(r1(t),meanlog=M, sdlog=S)) +log(dlnorm(r2(t),meanlog=M, sdlog=S)) +log(dlnorm(r3,meanlog=M, sdlog=S))+log(dlnorm(r4(t),meanlog=M, sdlog=S)) ;
 }
 
 densityFromRates <- function (t) {
   exp(logRateDensity(t)); 
 }
 
 
 A=integrate(densityFromRates,0,T)
 B=A$value
 
 G<- function (t) {
   densityFromRates(t)/B; 
 }