args = commandArgs(trailingOnly=TRUE)
source('EfficiencyCompare_utils.R')
#source('~/Desktop/validation/efficiency/simulated/EfficiencyCompare_utils.R')

n.sim <- args[1]
data.name <- args[2]
folder.path <- args[3]


#n.sim <- 20
#data.name <- "RSV2"
#folder.path <- "~/Desktop/validation/efficiency/others/"
  
Model <- c("Cons","Category")

for (model in Model) {
  
simulated.folder = paste0(folder.path,data.name,model,"/")

if (data.name=="RSV2" || data.name=="Shankarappa" || data.name == "anolis") {
  if (model == "Cons") {
    model.Nr = 2
  } else {
    model.Nr = 1
  }
  output.txt.folder = paste0(simulated.folder,"output/output_",data.name, "_", model.Nr,"_")
} else {
  output.txt.folder = paste0(simulated.folder,"output/output_",data.name, model,"20taxa_")
}

#calculation time in screen log file
Time.df <- c()
for (sim in 1:n.sim) {
  output.txt = readLines(paste0(output.txt.folder, sim, ".txt"))
  timeLine = output.txt[grep("Total calculation time", output.txt)]
  Time.df[sim] = as.numeric(gsub(" .+", "", gsub(".+[:] ", "", timeLine))) / 60 
}

#read loganalyser output of all simulations
ESS.txt <- read.table(paste0(folder.path,"ess/ESS_",data.name,model,".txt"),sep="\t", header=T)

#ESS of paramaters of interest
if (data.name == "Short" || data.name == "Medium") {
assign(paste0(model,".Efficiency.df"), get.simulated.efficiency(ESS.txt,Time.df))
}

if (data.name == "RSV2") {
  assign(paste0(model,".Efficiency.df"), get.RSV2.efficiency(ESS.txt,Time.df))
}
if (data.name == "Shankarappa") {
  assign(paste0(model,".Efficiency.df"), get.Shankarappa.efficiency(ESS.txt,Time.df))
}
if (data.name == "anolis") {
  assign(paste0(model,".Efficiency.df"), get.anolis.efficiency(ESS.txt,Time.df))
}

}

ratio.df <- Cons.Efficiency.df/Category.Efficiency.df
Mean <- apply(ratio.df,2,mean,trim=0.05)
write.table(t(Mean),file=paste0(folder.path, "EfficiencyTable_",data.name,".txt"),quote=F,sep="\t",row.names = FALSE)

pdf (file=paste0(folder.path, "EfficiencyCompare_",data.name,".pdf"),width=7,height=4)
boxplot(Cons.Efficiency.df/Category.Efficiency.df, xlim = c(0,14), boxwex = 0.2, at = seq(0,14,length.out=14),outline=FALSE,xaxt="n",yaxt="n")
axis(1, at = seq(0,14,length.out=14), labels = colnames(Cons.Efficiency.df),las=2,cex.axis=0.8)
axis(2,cex.axis=0.8)
abline(h=1.0,col="red",lwd=2.0)
title(main=data.name, ylab = c("ESS(cons/categories)"))
Mean <- apply(ratio.df,2,mean,trim=0.05)
dev.off()
