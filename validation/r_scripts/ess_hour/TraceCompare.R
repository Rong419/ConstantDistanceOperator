#step1: for rates
# read data
P <- read.table(file="/Users/rzha419/Desktop/4.txt", sep="\t", header=T)

# get ESS
U=P[10,]

E = matrix(nrow = 5, ncol = 1);
for (j in 1:5){
  E[j]=U[j+1];
}

E4=matrix(nrow = 5, ncol = 1);
for (j in 1:5){
    X=factor(E[[j]]);
    E4[j]=levels(X);
}
par(new=TRUE)
plot(E4,xlab = "",type='o', ylab = "",xlim=c(1,6),ylim=c(50,550),col="green")

plot(E1,xlab = "",type='o', ylab = "",ylim=c(100,600),col="black")

#step2:for tree height
Q <- read.table(file="/Users/rzha419/Desktop/6.txt", sep="\t", header=T)

# get ESS
V=Q[10,]

F = matrix(nrow = 6, ncol = 1);
for (j in 1:6){
  F[j]=V[j+1];
}

F6=matrix(nrow = 6, ncol = 1);
for (j in 1:6){
  X=factor(F[[j]]);
  F6[j]=levels(X);
}
par(new=TRUE)
plot(F8,xlab = "",type='o', ylab = "",xlim=c(1,6),ylim=c(450,1450),col="green")

plot(F1,xlab = "",type='o', ylab = "",ylim=c(100,600),col="black")