##This R code shows the failure of convergence 
##if the Nadaraya-Watson regression estimator and its gradient
##are directly used with a mean shift idea
##This simulation is used to justify the statement in Remark 2.1.
#############

##The function msOperatorNW() returns one iteration of the mean shift using the Nadaraya-Watson regression estimator  
## The Epanechnikov kernel is used as the kernel g.
##Inputs: 
##x_0: current location; 
##x: covariate observations
##y: response observations
##bw: bandwidth
msOperatorNW <- function(x_0,x,y,bw){
  radi <- bw
  ne.id <- ball.id(x,x_0,radi)
  res <- c()
  if (length(ne.id)>0){
    ne <- x[ne.id,,drop=FALSE]
    n.ne <- length(ne.id)
    dif <- sweep(as.matrix(ne),2,as.matrix(x_0))
    dis <- sqrt(rowSums(dif^2))
    weig <- 1-(dis/bw)^2
    weik <- (1-(dis/bw)^2)^2
    
    weit <- y[ne.id]*weig*sum(weik) - weig*sum(y[ne.id]*weik)
    n <- nrow(x)
    res <- colSums(weit*x[ne.id,,drop=FALSE])/sum(weit)
  }else{
    res <- c(NA,NA)
  }
  res
}

##msAlgorithmUpdate() implements one iteration of the regression mean shift algorithmtheusing the Nadaraya-Watson regression estimator at multiple locations x0
##Inputs: 
##x0: current locations; 
##x: covariate observations
##y: response observations
##bw: bandwidth
msAlgorithmUpdateNW <- function(x0,x,y,bw,kernel){
  M <- x0
  n <- nrow( x0 )
  
  for( i in 1:n ){
    if (is.na(x0[i,])){
      M[i,] <- x0[i,]
    }else{
      M[i,] <- msOperatorNW(x_0=x0[i,], x=x,y=y,bw=bw)
    }
  }
  output <- M
  return( output )
}

#############
library(mvtnorm)
library(tmvtnorm)
library(np)
library(rgl)

source("RMS_funs.r")

set.seed(1234)

n <- 200
bw <- 1.5
trans_method <- 1

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


xlim=seq(-3,3,length.out = 200)
ylim=seq(-3,3,length.out = 200)

fout <- outer(xlim, ylim, Vectorize(func))

image(xlim,ylim,fout, col=terrain.colors(40), ann=FALSE)
points(x)

modes.old <- x

tol.stop=1e-6
stop.id <- rep(0,nrow(modes.old))
k <- 1

modes.final <- modes.old
track.id <- 1:nrow(modes.old)
while(sum(!stop.id)>0){
  cat(k,",",sum(!stop.id),"\n")
  modes.new <- msAlgorithmUpdateNW(modes.old,x,yt,bw)##Regression Mean Shift using a bandwidth bw
  segments(modes.old[,1],modes.old[,2],modes.new[,1],modes.new[,2],col="red")
  stop.id <- (sqrt(rowSums((modes.old - modes.new)^2))<tol.stop)
  stop.id <- stop.id | is.na(stop.id)
  modes.old <- modes.new[!stop.id,,drop=FALSE]
  if (sum(stop.id)>0){
    delete.id <- track.id[stop.id]
    track.id <- track.id[!stop.id]
    modes.final[delete.id,] <- modes.new[stop.id,]
  }
  k <- k+1
}
##There is an error because the mean shift based on Nadaraya-Watson regression estimator
##easily moves points to regions where no data points are observed within a radius of the bandwidth h.

points(modes.final[,1],modes.final[,2],col="blue",pch=4)

tol.epsilon=1e-3

output <- connectedComponents( X=t(modes.final),
                               tol.epsilon=tol.epsilon )

plot(data.full, col=output$labels+2, cex=0.8,
     pch=16, main="h=1.5", xlab="x", ylab="y")
points(modes.final,pch=4)



