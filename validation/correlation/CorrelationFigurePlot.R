library(corrgram)
library(ape)

args = commandArgs(trailingOnly=TRUE)

log.file.path <- args[1]
tree.file.folder <- args[2]
output.figure.folder <- args[3]

#log.file.path <- "~/Desktop/validation/correlation/"
#tree.file.path <- "~/Desktop/validation/correlation/ratites/"
#output.figure.folder <- "~/Desktop/validation/correlation/"

ExtRates.log <- read.table(file = paste0(log.file.path, "ExtRates.log"), sep = "\t", header = T)
IntRates.log <- read.table(file = paste0(log.file.path, "IntRates.log"), sep = "\t", header = T)

External.Rates <- ExtRates.log[,2:8]
Internal.Rates <- IntRates.log[,2:6]
TreeHeight <- ExtRates.log$Tree.Height
t.AD <-  ExtRates.log$tMRCA.AD.
t.CE <- ExtRates.log$tMRCA.CE.
t.CEK <- ExtRates.log$tMRCA.CEK.
t.CEKO <- ExtRates.log$tMRCA.CEKO.
t.ADKCEKO <- ExtRates.log$tMRCA.ADCEKO.
Monophyly <- ExtRates.log[,9:13]

data.df <- cbind(External.Rates, Internal.Rates, t.AD, t.CE, t.CEK, t.CEKO, t.ADKCEKO, TreeHeight, Monophyly)
S1 <- data.df[which(data.df$Monophyly.AD.==1) ,]
S2 <- S1[which(S1$Monophyly.ADCEKO.==1),]
S3 <- S2[which(S2$Monophyly.CE.==1),]
S4 <- S3[which(S3$Monophyly.CEK.==1),]
S5 <- S4[which(S4$Monophyly.CEKO.==1),]

plot.df <- S5[,-(19:23)]
names(plot.df)<-c("r1","r2","r3","r4","r5","r6","r7","r8","r9","r10","r11","r12","t1","t2","t3","t4","t5","T")

pdf(paste(output.figure.folder,"RatesAndTimeCorrelation.pdf"), height=10, width=10)
corrgram(plot.df)
dev.off()

s.trees <- read.nexus(paste0(tree.file.path,"ratites.subst.trees"))

d1 <- sapply(strees, function(x) {x$edge.length[1]})
d2 <- sapply(strees, function(x) {x$edge.length[2]})
d3 <- sapply(strees, function(x) {x$edge.length[3]})
d4 <- sapply(strees, function(x) {x$edge.length[4]})
d5 <- sapply(strees, function(x) {x$edge.length[5]})
d6 <- sapply(strees, function(x) {x$edge.length[6]})
d7 <- sapply(strees, function(x) {x$edge.length[7]})
d8 <- sapply(strees, function(x) {x$edge.length[8]})
d9 <- sapply(strees, function(x) {x$edge.length[9]})
d10 <- sapply(strees, function(x) {x$edge.length[10]})
d11 <- sapply(strees, function(x) {x$edge.length[11]})
d12 <- sapply(strees, function(x) {x$edge.length[12]})

distance.df <- cbind(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12)
pdf(paste(output.figure.folder,"DistanceCorrelation.pdf"), height=10, width=10)
corrgram(distance.df)
dev.off()
