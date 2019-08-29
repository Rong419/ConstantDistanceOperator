#This script is to validate the operations of SimpleDistance
#that proposes a root time 
#and maintains genetic distances of the two branches linked to the root
##############


d1 <- 0.1; d2 <- 0.2; d3 <- 0.4; d4 <- 0.27
t <- 1;

r1 <- 0.1;
r2 <- 0.2;
r3 <- function (T) { d3 / T;}
r4 <- function (T) { d4 / (T-t);}

M = -3;
S = 0.25;

logRateDensity <- function (T) {
logd <- log(dlnorm(r1,meanlog=M, sdlog=S)) + log(dlnorm(r2,meanlog=M, sdlog=S)) + log(dlnorm(r3(T),meanlog=M, sdlog=S))+log(dlnorm(r4(T),meanlog=M, sdlog=S)) ;
}

densityFromRates <- function (T) {
exp(logRateDensity(T)); 
}

A=integrate(densityFromRates,1,15)
B=A$value

G<- function (T) {
 	densityFromRates(T)/B; 
 }
integrate(G,1,15)

E<- function (T) {
	 T*G(T); 
 }
 
integrate(E,1,15)
sqrt()

#minus the mean manually
T1<- function (T) {
  (T-7.818691)*(T-7.818691); 
 }
 
F<- function (T) {
  T1(T)*G(T); 
  }

 integrate(F,1,15)


l <- read.table(file="/Users/rzha419/Documents/SoftWare/beast2-master/dna.log", sep="\t", header=T)
 hist(l$mrcatime.ingroup.,prob=T,breaks=80)
  curve(G,from=0,to=T, col="red",add=TRUE)