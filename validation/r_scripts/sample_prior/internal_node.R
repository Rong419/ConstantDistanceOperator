libarary(ape)

args = commandArgs(trailingOnly=TRUE)


newick.tree.path <- args[1]
rates.path <- args[2]
logfile.folder <- args[3]
n.sim <- as.numeric(args[3])

newick.tree.path <- "/Users/rzha419/Workspace/ConstantDistanceOperator/validation/sample_prior/internal_nodes/test_internalnode_trees.txt"
rates.path <- "/Users/rzha419/Workspace/ConstantDistanceOperator/validation/sample_prior/internal_nodes/test_internalnode_rates.txt"
logfile.folder <- "/Users/rzha419/Workspace/ConstantDistanceOperator/out/artifacts/ConstantDistanceOperator_jar/" 

trees = readLines(newick.tree.path)
rates = read.table(file=rates.path)

Mean = c()
Stdev = c()

for (idx in 1:length(trees)){
tree = read.tree(text = paste0(trees[1],";"))
time = tree$edge.length
rate = as.numeric(rates[idx,])
#initialize the constant distance and root time
d1 = rate[1] * time[4]; 
d2 = rate[2] * time[3]; 
d3 = rate[3] * time[1]; 
d4 = rate[4] * time[2]
T = time[1];

#functions of rates on the 4 branches
#i.e. distance divided by the time duration
r1 <- function (t) { d1 / t;}
r2 <- function (t) {d2 / t;}
r3 = rate[3]
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
#make sure that N should be equal to 1
G <- function (t) {
 	densityFromRates(t)/B; 
 }
N = integrate(G,0,T)


#the mean of the distribution
E <- function (t) {
	 t*G(t); 
 } 
mean = integrate(E,0,T)


#the standard deviation of the distribution
Mean[idx] = mean$value
F <- function (t) {
  (t - Mean[idx]) * (t - Mean[idx]) * G(t); 
  }
  
variance = integrate(F, 0, T)
Stdev[idx] = sqrt(variance$value)
}


for (scenario in 1:length(trees)) {
     for (sim in 1:n.sim) {
#read log file
l <- read.table(file=paste0(logfile.folder,"internalnode_S",[ScenarioHere],"_",,".log"), sep="\t", header=T)

#plot the histogram of the parameter
#i.e. the distribution of sampled root time 
 hist(l$mrcatime.ingroup.,prob=T,breaks=80,xlab="",main='')

#add the curve of the theoretical distribution of the root time
#the normalized CDF distribution
curve(G,from=0,to=T, col="red",add=TRUE,lwd=2)
curve(densityFromRates,from=0,to=T, col="red",add=TRUE)
}
}
d <- density(l$mrcatime.ingroup., n=length(l$mrcatime.ingroup.))
plot(d, main="")

