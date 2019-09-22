library(ape)

args = commandArgs(trailingOnly=TRUE)


newick.tree.path <- args[1]
rates.path <- args[2]
logfile.folder <- args[3]
n.sim <- as.numeric(args[4])
output.figure.folder <- args[5]

#newick.tree.path <- "/Users/rzha419/Workspace/ConstantDistanceOperator/validation/sample_prior/internal_nodes/test_internalnode_trees.txt"
#rates.path <- "/Users/rzha419/Workspace/ConstantDistanceOperator/validation/sample_prior/internal_nodes/test_internalnode_rates.txt"
#logfile.folder <- "/Users/rzha419/Workspace/ConstantDistanceOperator/out/artifacts/ConstantDistanceOperator_jar/" 
#output.figure.folder <- "/Users/rzha419/Workspace/ConstantDistanceOperator/out/artifacts/ConstantDistanceOperator_jar/" 

trees = readLines(newick.tree.path)
rates = read.table(file=rates.path)

Mean = c()
Stdev = c()
Abs = c()
Table = c()

for (idx in 1:length(trees)){
tree = read.tree(text = paste0(trees[idx],";"))
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

#population size
N = 0.3

#distribution of rates
densityFromRates <- function (t) {
dlnorm(r1(t),meanlog=M, sdlog=S) * dlnorm(r2(t),meanlog=M, sdlog=S) * dlnorm(r3,meanlog=M, sdlog=S) * dlnorm(r4(t),meanlog=M, sdlog=S);
}

#distribution of times
densityFromTimes <- function (t) {
(1 / N) * exp(-3 * t / N) * (1 / N) * exp(-(T-t) / N);
}

#distribution of tree
densityFromTree <- function (t) {
densityFromTimes(t) * densityFromRates(t);
}

#integrate the density function
#the area of the curve
A=integrate(densityFromTree,0,T)
B=A$value

#normalize the density function
#make sure that Normalized should be equal to 1
G <- function (t) {
 	densityFromTree(t)/B; 
 }
Normalized = integrate(G,0,T)


#the mean of the distribution
E <- function (t) {
	 t*G(t); 
 } 
mean = integrate(E,0,T)
Abs[idx] = mean$abs.error

#the standard deviation of the distribution
Mean[idx] = mean$value
F <- function (t) {
  (t - Mean[idx]) * (t - Mean[idx]) * G(t); 
  }
  
variance = integrate(F, 0, T)
Stdev[idx] = sqrt(variance$value)



     Table[idx] = paste0("Scenario", idx," : M = ",Mean[idx], " with absolute error < ", Abs[idx], "S = ", Stdev[idx])
     
     
     for (sim in 1:n.sim) {
         pdf (paste0(output.figure.folder,"internlanode_S",idx,"_",sim,".pdf"),width=7,height=5)
         
         #read log file
         l <- read.table(file=paste0(logfile.folder,"internalnode_S",idx,"_",sim,".log"), sep="\t", header=T)

         #plot the histogram of the parameter
         #i.e. the distribution of sampled root time 
         hist(l$mrca.age.ingroup.,prob=T,breaks=80,xlab="",main='')

         #add the curve of the theoretical distribution of the root time
         #the normalized CDF distribution
         curve(G,from=0,to=T, col="red",add=TRUE,lwd=2)
         
         dev.off()
         }
                  
}
      write.table(Table,file=paste0(output.figure.folder,"numerical_internalnode.txt"))


