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

logRateDensity <- function (d3) {
logd <- log(dlnorm(r1,meanlog=M, sdlog=S)) + log(dlnorm(r2,meanlog=M, sdlog=S)) + log(dlnorm(r3(d3),meanlog=M, sdlog=S))+log(dlnorm(r4(d3),meanlog=M, sdlog=S)) ;
}

densityFromRates <- function (d3) {
exp(logRateDensity(d3)); 
}

A=integrate(densityFromRates,0,0.67)
B=A$value

G<- function (d3) {
 	densityFromRates(d3)/B; 
 }
integrate(G,0,0.67)

E<- function (d3) {
	 d3*G(d3); 
 }
 
integrate(E,0,0.67)


T1<- function (d3) {
  (d3-0.3476093)*(d3-0.3476093); 
 }
 
F<- function (d3) {
  T1(d3)*G(d3); 
  }

integrate(F,0,0.67)
sqrt(V)

#read tree files
strees <- read.nexus("/Users/rzha419/Documents/SoftWare/beast2-master/dna1.subst.trees")
#get the distance on branch
e <- sapply(strees, function(x) {x$edge.length[4]})
#plot the distribution of the distance
hist(e,prob=T, breaks=100)
#plot the integrated distribution
curve(G,from=0,to=0.67, col="red",add=TRUE)
curve(densityFromRates,from=0,to=0.67, col="red")