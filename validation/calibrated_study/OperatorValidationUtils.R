library(HDInterval)
library(ggplot2)
library(gtools)
library(ggpubr)
library(ape)
library(rowr)
library(grid)

get.95 <- function(a.vector) {
  res = hdi(a.vector, cresMass=.95)
  return(c(res[[1]], res[[2]], mean(a.vector)))
}

get.calibrated.plot <- function(true.df,log.df, n.sim, x.min, x.max, y.min, y.max,prior.mean) {
     full.df <- data.frame(cbind(true.df,log.df))
     colnames(full.df) <- c("true","lower","upper","mean")
     full.df[,1] <- jitter(full.df[,1],factor=1)
  
     plot.hdi <- matrix()
     T = 0
     F = 0
     for (i in 1:n.sim) {
        if (full.df[i,3] >= true.df[i] & true.df[i] >= full.df[i,2]) {
            plot.hdi[i] <- TRUE
            T <- T + 1
        } else {
           plot.hdi[i] <- FALSE
           F <- F + 1
        }
     }
  
     reg.df = data.frame(cbind(full.df,plot.hdi))
     false.df <- reg.df[!plot.hdi,]
     true.df <- reg.df[plot.hdi,]
     
     p <- ggplot() + 
     geom_linerange(false.df,mapping=aes(x=false.df$true, ymax=false.df$upper, ymin=false.df$lower), color="red", alpha=.4, size=1) +
     geom_linerange(data=true.df[plot.hdi,],mapping=aes(x=true.df$true, ymax=true.df$upper, ymin=true.df$lower), color="lightblue", alpha=.4, size=1) +
     geom_point(full.df,mapping=aes(x=full.df$true, y=full.df$mean),shape=20,size=2) + 
     geom_abline(slope=1, linetype="dotted") + 
     coord_cartesian(xlim=c(x.min, x.max),ylim=c(y.min, y.max)) +
     xlab("True value") + ylab("Posterior mean") +
     geom_abline(slope=0, intercept=prior.mean, color="blue") +
     geom_abline(slope=1, linetype="dotted") +
     scale_x_continuous(breaks=seq(x.min,x.max,1)) +
     scale_y_continuous(breaks=seq(y.min,y.max,1)) +
     theme(
       panel.grid.minor = element_blank(),
       panel.border = element_blank(),
       panel.background = element_blank(),
       plot.background = element_blank(),
       plot.title = element_text(hjust=0.5),
       axis.line = element_line(),
       axis.ticks = element_line(color="black"),
       axis.text.x = element_text(color="black", size=10),
       axis.text.y = element_text(color="black", size=10),
       axis.title.x = element_text(size=12),
       axis.title.y = element_text(size=12)
     )
   return (p)
}

