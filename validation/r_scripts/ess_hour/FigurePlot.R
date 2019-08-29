#This R script is used to compare the ESS of the parameters
#Series 5
#To plot the necessary figures that show the comparisons

#show the ESS of parameters in "M1 M2 M3 M4" matrix
plot(M1,type='o',xlab = "parameter",ylim=c(0,3400), ylab = "ESS",xaxt="n")
par(new=TRUE)
plot(M2,xlab = "",type='o', ylab = "",ylim=c(0,3400), col="red",xaxt="n")
par(new=TRUE)
plot(M3,xlab = "",type='o', ylab = "",ylim=c(0,3400), col="blue",xaxt="n")
par(new=TRUE)
plot(M4,xlab = "",type='o', ylab = "",ylim=c(0,3400), col="green",xaxt="n")
#show the name of parameters on x-axim 
axis(1,at=seq(1,14,1),labels=c("posterior","likelihood","prior","height","r_m","r_v","r_c","kapa","birth","death","q1","q2","q3","q4"))


#show comparison by X vs Y
plot(M1,M3,xlim=c(0,2000),ylim=c(0,2000),col="red",xlab="withcons",ylab="without cons/quantiles",xaxt="n")
par(new=TRUE)
plot(M1,M4,xlim=c(0,2000),ylim=c(0,2000),col="blue",xlab="",ylab="")

#show the equal line, X=Y
criterion <- function(m) {
  m;
} 
curve(criterion,0, 2000, add = TRUE, col = "black",lwd=1,lty=2)

#show the fitted line of X vs Y
fit1<-lm(M3~M1)
abline(fit1,col = "red", xlim=c(0,2000), ylim=c(0,2000),lwd=2)

fit2<-lm(M4~M1)
abline(fit2,col = "blue", xlim=c(0,2500), ylim=c(0,2500),lwd=2)


#construct data matrix
d <- data.frame(po=posterior1,Li=likelihood1,pr=prior1,tr=treelikelihood1,he=treeheight1,rm=ratemean1,rv=ratevariance1,rc=ratecoef1,k=kapa1,b=birthrate1,de=deathrate1)
f <- data.frame(po=posterior2,Li=likelihood2,pr=prior2,tr=treelikelihood2,he=treeheight2,rm=ratemean2,rv=ratevariance2,rc=ratecoef2,k=kapa2,b=birthrate2,de=deathrate2)
g <- data.frame(po=posterior3,Li=likelihood3,pr=prior3,tr=treelikelihood3,he=treeheight3,rm=ratemean3,rv=ratevariance3,rc=ratecoef3,k=kapa3,b=birthrate3,de=deathrate3)

#show box plot
boxplot(posterior1,likelihood1,prior1,treeheight1,ratemean1,kapa1,birthrate1,deathrate1,freq11,freq12,freq13,freq14,xlim=c(1,12),ylim=c(0,3000),outline=FALSE,boxwex=0.2,at=1:12 - 0.3)
boxplot(posterior2,likelihood2,prior2,treeheight2,ratemean2,kapa2,birthrate2,deathrate2,freq21,freq22,freq23,freq24,col="green",add=TRUE,xlim=c(1,12),ylim=c(0,3000),outline=FALSE,boxwex=0.2,at=1:12)
boxplot(posterior3,likelihood3,prior3,treeheight3,ratemean3,kapa3,birthrate3,deathrate3,freq31,freq32,freq33,freq34,col="blue",add=TRUE,xlim=c(1,12),ylim=c(0,3000),outline=FALSE,boxwex=0.2,at=1:12+0.3)
abline(v=1.4,col="green")

#add legend
legend("top",inset=.02,c("categories"),pch=2,col="black",box.lty = 0)
legend("top",inset=.08,c("cons"),pch=1,col="red",box.lty = 0)
legend("top",inset=.15,c("nocons"),pch=16,col="green",box.lty = 0)
legend("topleft",inset=.01,c("cons","categories","nocons"),pch=1,col=c(1,2,4),box.lty = 0)
