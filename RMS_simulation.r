library(mvtnorm)
library(tmvtnorm)
library(np)
library(rgl)

source("RMS_funs.r")

for (sim_num in 1:200){
  for (n in c(200,500,1000)){
    for (trans_method in c(1,2)){
      
      set.seed(762*sim_num + 121231)
      
      i.sample <- sim_num
      
      trans <- function(x) 1/(1+exp(-10*x)) + 0.01
      
      func <- function(x,y) ##regression function in the model
        dmvnorm(cbind(x,y),c(1,1),diag(0.5,2)) + dmvnorm(cbind(x,y),c(-1,-1),diag(c(0.3,0.9))) 
      
      x <- rtmvnorm(n,c(0,0),sigma=diag(1.5,2),lower=c(-2,-2),upper=c(2,2))
      y <- func(x[,1],x[,2])  + rnorm(n,0,0.1)
      
      if(trans_method==1){##Transformation 1
        yt <- trans(y) 
      }else if(trans_method==2){##Transformation 2
        if (min(y)<0.1){
          yt <- y - min(y) + 0.1
        }
      }
      
      bw <- gbbw(x,yt)##bandwidth selection
      bw 
      
      modes.final <- RegMS(x,yt,bw)##Regression Mean Shift using the selected bandwidth
      
      output <- connectedComponents( X=t(modes.final))
      
      modes.est <- t(output$components)
      mode.num <- NROW(modes.est)
      
      
      bwseq <- seq(0.9,2.5,by=0.1)##Regression Mean Shift using bandwidths on a grid
      multimodality <- c()
      for (bw in bwseq){
        modes.seq <- RegMS(x,yt,bw)
        
        output.seq <- connectedComponents( X=t(modes.seq))
        
        modes.est.seq <- t(output.seq$components)
        multimodality <- c(multimodality,NROW(modes.est.seq))
      }
      
      result <- list(mode.num=mode.num,modes.est=modes.est,multimodality=multimodality)
      
      save(result,file=paste("/output_",i.sample,".RData",sep=""))
    }
  }
}
