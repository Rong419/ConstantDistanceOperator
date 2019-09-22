args = commandArgs(trailingOnly=TRUE)

source('OperatorValidationUtils.R')

n.sim <- as.numeric(args[1])
taxa <- args[2]
log.file.path <- args[3]
true.txt.path <- args[4]
output.figure.folder <- args[5]

#log.file.path <- "~/Desktop/validation/calibrated/logs/"
#true.txt.path <- "~/Desktop/validation/calibrated/true/"
#output.figure.folder <- "~/Desktop/validation/calibrated/figures/"

# the prefix of true data file .txt
# which is also the abbreviations of parameters to compare
parameter.names <- c("TreH","Kap","Ucld","Bir")

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
   log.Bir[i,1:3] = get.95(this.sim.df$BirthRate[501:5001])
}

TreH.figure <- get.calibrated.plot(true.TreH,log.TreH, n.sim, 0,4,0,4,"TreeHeight")
kap.figure <- get.calibrated.plot(true.Kap,log.Kap, n.sim, 0,1,0,1,"Kappa")
Ucld.figure <- get.calibrated.plot(true.Ucld,log.Ucld, n.sim, 0,0.6,0,0.6,"UcldStdev")
BirthRate.figure <- get.calibrated.plot(true.Bir,log.Bir, n.sim,1,4,1,4,"BithRate")

Figure <- ggarrange(TreH.figure,kap.figure ,Ucld.figure ,BirthRate.figure ,ncol=2,nrow=2)

pdf(paste(output.figure.folder,"CalibratedCoverage.pdf"), height=8, width=8)
annotate_figure(Figure,top = text_grob("Compare of true and estimated values in well-calibrated simulation study",face = "bold", size = 14))
dev.off()


