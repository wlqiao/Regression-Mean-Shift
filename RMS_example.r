##This R code reproduces Figure 4.2

source("RMS_funs.r")
#############
library(mvtnorm)
library(tmvtnorm)
library(np)
library(rgl)
set.seed(12345)

n <- 200
trans_method <- 1

trans <- function(x) 1/(1+exp(-10*x)) + 0.01

func <- function(x,y)
  dmvnorm(cbind(x,y),c(1,1),diag(0.5,2)) + dmvnorm(cbind(x,y),c(-1,-1),diag(c(0.3,0.9))) 

x <- rtmvnorm(n,c(0,0),sigma=diag(1.5,2),lower=c(-2,-2),upper=c(2,2))
y <- func(x[,1],x[,2])  + rnorm(n,0,0.1)

if(trans_method==1){
  yt <- trans(y)
}else if(trans_method==2){
  if (min(y)<0.1){
    yt <- y - min(y) + 0.1
  }
}

bw <- gbbw(x,yt)
bw 

##########################


par(mfrow=c(2,2))
modes.final <- RegMS(x,yt,bw,show.graph=TRUE,background=TRUE,x1lim=c(-2.3,2.3),x2lim=c(-2.3,2.3))

bw_test <- 1
modes <- RegMS(x,yt,bw_test,x1lim=c(-2,-2),x2lim=c(-2,-2))
output <- connectedComponents( X=t(modes))
col.names <- c("green","blue","pink","cyan","brown","orange","grey50","red","lightcoral")

plot(x, col=col.names[output$labels], cex=0.8,
     pch=16, xlab="(b)", ylab="", main=paste("h=",bw_test,sep=""))
points(modes,pch=4)


modes <- RegMS(x,yt,bw,x1lim=c(-2.3,2.3),x2lim=c(-2.3,2.3))
output <- connectedComponents( X=t(modes))
col.names <- c("green","blue")

plot(x, col=col.names[output$labels], cex=0.8,
     pch=16, xlab="(c)", ylab="",main=paste("h=",round(bw,1),sep=""))
points(modes,pch=4)

bw_test <- 2.5
modes <- RegMS(x,yt,bw_test,x1lim=c(-2.3,2.3),x2lim=c(-2.3,2.3))
output <- connectedComponents( X=t(modes))
col.names <- c("green")

plot(x, col=output$labels+2, cex=0.8,
     pch=16, main=paste("h=",bw_test,sep=""),xlab="(d)")
points(modes,pch=4)
