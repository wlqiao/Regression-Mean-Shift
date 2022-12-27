##Aggregate the simulation results from the output of RMS_simuation.r
##This R code reproduces Figure 4.1

library(pracma)
n <- 200
#n <- 500
#n <- 1000
par(mfrow=c(1,2))

for (method in c(1,2)){
  cat("method ", method,"\n")
  mode.num <- c()
  multimodality <- c()
  for (i in 1:200){
    try(load(file=paste("T_",method,"_n_",n,"/output_",i,".RData",sep="")))
    mode.num <- c(mode.num, result$mode.num)
    multimodality <- rbind(multimodality, result$multimodality)
  }
  
  bimode_rate <- length(which(mode.num==2))/200
  cat(bimode_rate,"\n")
  
  colnames(multimodality)=seq(0.9,2.5,by=0.1)
  
  perc <- apply(multimodality,2, function(x) length(which(x==2)))/200
  
  par(mar = c(5,5,2,5))
  boxplot.matrix(multimodality,xlim=c(1,17),ylim=c(0,15),xlab="h",ylab="number of modes",main=paste("Transformation T",method,sep="")) 
  
  library(scales)
  
  par(new = T)
  plot(perc, ylim=c(0,0.9), type="o", lty=2, axes=F, xlab=NA, ylab=NA, col="red",xlim=c(1,17))
  yticks_val <- pretty_breaks(n=5)(perc)
  axis(side = 4,at=yticks_val, lab=percent(yticks_val))
  mtext(side = 4, line = 3, 'relative frequency of bimodality')
  
}
