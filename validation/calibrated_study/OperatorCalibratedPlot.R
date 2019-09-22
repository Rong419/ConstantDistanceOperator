args = commandArgs(trailingOnly=TRUE)

source('OperatorValidationUtils.R')

n.sim <- args[1]
taxa <- args[2]
log.file.path <- args[3]
true.txt.path <- args[4]
# the prefix of true data file .txt
# which is also the abbreviations of parameters to compare
parameter.names <- c("TreH","Kap","Ucld")

# read true data from .txt files
for (para.name in parameter.names) {
     assign(paste0("true.",para.name),read.table(file=paste0(true.txt.path,taxa,"taxa/",para.name,".txt"), sep="\t"))
     
     # intilaize the data frame to save the sampled data from .log files
     assign(paste0("log.",para.name), data.frame(matrix(ncol=3, nrow=n.sim)))
}

# iterating all the simulations to get 95% HPD and posterior mean
for (i in 1:n.sim) {
   this.sim.df = read.table(file=paste0(log.file.path,"/calibrated",taxa,"taxa_",i,".log"),sep="\t",header=T)
   log.TreH[i,1:3] = get.95(this.sim.df$Tree.height[501:5001])
   log.Kap[i,1:3] = get.95(this.sim.df$kappa[501:5001])
   log.Ucld[i,1:3] = get.95(this.sim.df$ucldStdev[501:5001])
}

TreH.figure <- get.calibrated.plot(true.TreH,log.TreH, n.sim, x.min, x.max, y.min, y.max,prior.mean)
kap.figure <- get.calibrated.plot(true.Kap,log.Kap, n.sim, x.min, x.max, y.min, y.max,prior.mean)
Ucld.figure <- get.calibrated.plot(true.Ucld,log.Ucld, n.sim, x.min, x.max, y.min, y.max,prior.mean)

