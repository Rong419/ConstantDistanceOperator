args = commandArgs(trailingOnly=TRUE)
source('EfficiencyCompare_utils.R')
#source('~/Desktop/validation/efficiency/EfficiencyCompare_utils.R')

n.sim <- args[1]
taxa <- args[1]
data.file.path <- args[3]
output.figure.folder <- args[4]

#n.sim <- 20
#taxa <- 20
#data.file.path <- "~/Desktop/validation/efficiency/simulated/"
  
n.model <- c("Cons","Category")
l.sequence <- c("Short","Medium")
n.taxa <- paste0(taxa,"taxa")


  for (sequence.length in l.sequence) {
    for (model.name in n.model) {
      output.txt.folder = paste0(simulated.folder,"output/output_", sequence.length, model.name, taxa.number)
      
      #calculation time in screen log file
      Time.df <- c()
      for (sim in 1:n.sim) {
      output.txt = readLines(paste0(output.txt.folder, sim, ".txt"))
      timeLine = output.txt[grep("Total calculation time", output.txt)]
      Time.df[sim] = as.numeric(gsub(" .+", "", gsub(".+[:] ", "", timeLine))) / 60 
      }

      #read loganalyser output of all simulations
      ESS.txt <- read.table(paste0(data.file.path,"ess/ESS_",sequence.length, model.name, taxa.number,".txt"),sep="\t", header=T)

      #ESS of paramaters of interest
      assign(paste0(model.name,".Efficiency.df"), get.simulated.efficiency(ESS.txt,Time.df))
    }
    assign(paste0(sequence.length,".Ratio.df"), Cons.Efficiency.df/Category.Efficiency.df)
    Mean <- apply(Cons.Efficiency.df/Category.Efficiency.df, 2, mean, trim=0.05)
    write.table(t(Mean),file=paste0(output.figure.folder, "EfficiencyTable_",data.name,".txt"),quote=F,sep="\t",row.names = FALSE)
  }
  #ShortRatio.df <- ShortCons.Efficiency.df/ShortCategory.Efficiency.df
  #MediumRatio.df <- MediumCons.Efficiency.df/MediumCategory.Efficiency.df


pdf (file=paste0(output.figure.folder, "EfficiencyCompare_",data.name,".pdf"),width=7,height=4)
boxplot(Cons.Efficiency.df/Category.Efficiency.df, xlim = c(0,14), boxwex = 0.2, at = seq(0,14,length.out=14),outline=FALSE,xaxt="n",yaxt="n")
axis(1, at = seq(0,14,length.out=14), labels = colnames(Cons.Efficiency.df),las=2,cex.axis=0.8)
axis(2,cex.axis=0.8)
abline(h=1.0,col="red",lwd=2.0)
title(main=data.name, ylab = c("ESS(cons/categories)"))
Mean <- apply(ratio.df,2,mean,trim=0.05)
dev.off()
