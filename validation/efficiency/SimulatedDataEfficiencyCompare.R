args = commandArgs(trailingOnly=TRUE)
source('EfficiencyCompare_utils.R')
#source('~/Desktop/validation/efficiency/EfficiencyCompare_utils.R')

n.sim <- args[1]
taxa <- args[2]
data.file.path <- args[3]
output.figure.folder <- args[4]

#n.sim <- 20
#taxa <- 20
#data.file.path <- "~/Desktop/validation/efficiency/simulated/"
#output.figure.folder <- "~/Desktop/validation/efficiency/figures/"
  
n.model <- c("Cons","Category")
l.sequence <- c("Short","Medium")
n.taxa <- paste0(taxa,"taxa")


  for (sequence.length in l.sequence) {
    
    for (model.name in n.model) {
      
      #calculation time in screen log file
      Time.df <- c()
      for (sim in 1:n.sim) {
      output.txt = readLines(paste0(data.file.path, sequence.length, model.name, "/output/output_", sequence.length, model.name, n.taxa, "_", sim, ".txt"))
      timeLine = output.txt[grep("Total calculation time", output.txt)]
      Time.df[sim] = as.numeric(gsub(" .+", "", gsub(".+[:] ", "", timeLine))) / 60 
      }

      #read loganalyser output of all simulations
      ESS.txt <- read.table(paste0(data.file.path,"ess/ESS_",sequence.length, model.name, n.taxa,".txt"),sep="\t", header=T)

      #ESS of paramaters of interest
      assign(paste0(model.name,".Efficiency.df"), get.simulated.efficiency(ESS.txt,Time.df))
    }
    assign(paste0(sequence.length,".Ratio.df"), Cons.Efficiency.df/Category.Efficiency.df)
    Mean <- apply(get(paste0(sequence.length,".Ratio.df")), 2, mean, trim=0.05)
    write.table(t(Mean),file=paste0(output.figure.folder, "EfficiencyTable_",sequence.length,".txt"),quote=F,sep="\t",row.names = FALSE)
  }
  #ShortRatio.df <- ShortCons.Efficiency.df/ShortCategory.Efficiency.df
  #MediumRatio.df <- MediumCons.Efficiency.df/MediumCategory.Efficiency.df


pdf (file=paste0(output.figure.folder, "EfficiencyCompare_Simulated",n.taxa,".pdf"),width=7,height=4)
boxplot(Short.Ratio.df, xlim = c(0,15), ylim=c(0,10), col="red", boxwex = 0.4, at = seq(0,15,length.out=14) - 0.2,outline=FALSE,xaxt="n",yaxt="n")
axis(1, at = seq(0,15,length.out=14), labels = colnames(Cons.Efficiency.df),las=2,cex.axis=0.8)
axis(2,cex.axis=0.8)
boxplot(Medium.Ratio.df, xlim = c(0,15), ylim=c(0,10),boxwex = 0.4, at = seq(0,15,length.out=14) + 0.2,col="green",outline=FALSE,xaxt="n",yaxt="n",add=TRUE)
abline(h=1.0,col="red",lwd=2.0)
legend("topright",inset=.01,c("Short","Medium"),fill =c("red","green"),box.lty = 0)
title(main=paste0("Simulated data set (",n.taxa, ")"), ylab = c("ESS(cons/categories)"))
dev.off()
