#This is an auxiliary script for comparing ESS
#by reading data from Excel file

#"ESS.xlsx" contains 3 columns for ESS obtained by different methods
library(readxl)
A <- read_excel("C:/Users/RongZhang/Desktop/ESS.xlsx");
#get values in ESS
ESS1 <- A$ESS1  #standard BEAST ussing categories
ESS2 <- A$ESS2  #standard BEAST ussing rates
ESS3 <- A$ESS3  #with ConstantDistance Operator
#another way
ESS1 = matrix(nrow = 37, ncol = 1);
for (i in 1:37){
  ESS1[i]=A[i+1];
}

#plot X vs Y 
plot(ESS1,ESS2,xlab = "", ylab = "", col = "red", xlim=c(0,8000), ylim=c(0,8000),pch=1)
par(new=TRUE)
plot(ESS1,ESS3,axes = FALSE,xlim=c(0,8000), ylim=c(0,8000),xlab = "", ylab = "", col = "green",pch=1)
par(new=TRUE)
plot(ESS2,ESS3,axes = FALSE,xlim=c(0,8000), ylim=c(0,8000),xlab = "", ylab = "", col = "blue",pch=1)

#fitted line 
fit1  <- lm(ESS1~ESS2)
fit2  <- lm(ESS1~ESS3)
fit3  <- lm(ESS2~ESS3)
abline(fit1,col = "red", xlim=c(0,8000), ylim=c(0,8000),lwd=2)
abline(fit2,col = "green",xlim=c(0,8000), ylim=c(0,8000),lwd=2)
abline(fit3,col = "blue",xlim=c(0,8000), ylim=c(0,8000),lwd=2)

legend("topleft", inset=.05, title="ESS Comparison",c("Ca vs Ra","Ca vs Dis","Ra vs Dis"), fill= rainbow(3), horiz=TRUE)

plot(ESS1,ESS2,
    main = "Result comparison of ESS in standard BEAST and ConstantDistance Operator",
     xlab = "standard BEAST",
     ylab = "ConstanttDistance Operator",
     col = "blue", 
    pch=1,
    axes = FALSE,
    xlim=c(0,15),
    ylim=c(0,15)
    )
par(new=TRUE)

#criterion line y = x
criterion <- function(m) {
  m;
} 
curve(criterion,0, 8000, add = TRUE, col = "black",lwd=1,lty=2)