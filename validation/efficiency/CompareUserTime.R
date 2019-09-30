args = commandArgs(trailingOnly=TRUE)

main.path <- args[1]
output.figure.folder  <- args[2]

# main.path <- "~/Desktop/efficiency/"
# output.figure.folder<- "~/Desktop/efficiency/figures"

data.name <- c("anolis", "RSV2")
# how many character to get in the .txt
time.length <- c(7, 8)

for (i in length(data.name)) {
     data = data.name[i]
     length = time.length[i]
     
     L =  readLines(paste0(main.path,data,"/time_",data,".txt"))

     categories = as.numeric(substr(L[1:100],1,length))
     
     cons = as.numeric(substr(L[101:200],1,length))
     
     # to show the mean time in the figure
     A = c(1.5,2.5)
     B = c(mean(categories),mean(cons))
     
     # set the y axis
     y.min = min(min(categories),min(cons))
     y.max = max(max(categories),max(cons))
     
     # make the figure
     pdf(paste0(output.figure.folder, "RunTimeCompare_", data, ".pdf"), width = 7, height = 5)
     boxplot(categories, ylim = c(y.min,y.max), xlim = c(1,3), boxwex = 0.8, at = 1.5)
     boxplot(cons, ylim = c(y.min, y.max), col="green", add = TRUE, xlim = c(1, 3), boxwex = 0.8, at = 2.5)
     axis(1, at = c(1.5, 2.5), labels = c("categories", "cons"))
     title(main = paste0("Runing time of ", data), ylab = "time (second)")
     par(new = TRUE)
     plot(A, B, ylim = c(5200,6500), xlim=c(1,3), type = 'o', xlab = '', ylab = '', col = "orange", xaxt = "n", yaxt = 'n')
     text(A, B, labels = B, pos= 4, cex = 0.8)
     dev.off()
}
