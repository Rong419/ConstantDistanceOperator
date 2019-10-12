args = commandArgs(trailingOnly=TRUE)
source('EfficiencyCompare_utils.R')
#source('~/Desktop/validation/efficiency/EfficiencyCompare_utils.R')

n.sim <- args[1]
data.name <- args[2]
data.file.path <- args[3]
output.figure.folder <- args[4]


#n.sim <- 100
#data.name <- "RSV2"
#data.file.path <- "~/Desktop/validation/efficiency/others/"
#output.figure.folder <- "~/Desktop/validation/efficiency/figures/"
  
Model <- c("Cons","Category")

for (model in Model) {

	simulated.folder = paste0(data.file.path,data.name,model,"/")

	if (model == "Cons") {
		model.Nr = 2
	} else {
		model.Nr = 1
	}
	output.txt.folder = paste0(simulated.folder,"output/output_",data.name, "_", model.Nr,"_")


	#calculation time in screen log file
	Time.df <- c()
	for (sim in 1:n.sim) {
		output.txt = readLines(paste0(output.txt.folder, sim, ".txt"))
		timeLine = output.txt[grep("Total calculation time", output.txt)]
		Time.df[sim] = as.numeric(gsub(" .+", "", gsub(".+[:] ", "", timeLine))) / 60 
	}

	#read loganalyser output of all simulations
	ESS.txt <- read.table(paste0(data.file.path,"ess/ESS_",data.name,model,".txt"),sep="\t", header=T)

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
write.table(t(Mean),file=paste0(output.figure.folder, "EfficiencyTable_",data.name,".txt"),quote=F,sep="\t",row.names = FALSE)

pdf (file=paste0(output.figure.folder, "EfficiencyCompare_",data.name,".pdf"),width=7,height=4)
boxplot(Cons.Efficiency.df/Category.Efficiency.df, xlim = c(0,14), boxwex = 0.5, at = seq(0,14,length.out=14),outline=FALSE,xaxt="n",yaxt="n")
axis(1, at = seq(0,14,length.out=14), labels = colnames(Cons.Efficiency.df),las=2,cex.axis=0.8)
axis(2,cex.axis=0.8)
abline(h=1.0,col="red",lwd=2.0)
title(main=data.name, ylab = c("Ratio of ESS per hour (cons/categories)"))
Mean <- apply(ratio.df,2,mean,trim=0.05)
dev.off()
