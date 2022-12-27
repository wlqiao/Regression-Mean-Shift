##ball.id() returns of the row numbers of data in the ball with center x and radius h
ball.id <- function(data,x,h){
  which(apply(data, 1, function(u) sum((u-x)^2))<=h^2)
}

##kde() returns the kernel density estimates at x using dat as sample and bandwidth h
##The kernel function used in kde() is Quartic.
kde <- function(dat,h,x){
  k <- nrow(x)
  n <- nrow(dat)
  res <- c()
  kern <- function(x) (abs(x)<=1)*3/pi*(1-x^2)^2
  
  for (i in 1:k){
    dif <- sweep(as.matrix(dat),2,as.matrix(x[i,,drop=FALSE]))
    dis <- sqrt(rowSums(dif^2))/h
    res <- c(res, 1/(n*h^2)*sum(kern(dis)))
  }
  res
}

##gbbw() returns the selected bandwidth for the regression mean shift algorithm.
##The bandwidth selection strategy is based on least square cross validation as described in Section 4 of the manuscript.
gbbw <- function(x,yt){
  
  n <- nrow(x)
  
  bw0 <- npregbw(formula = yt ~ x[,1] + x[,2], ckertype="epanechnikov")
  
  bw1 <- bw0
  bw1$bw <- bw0$bw * n^(1/48)
  
  kre1 <- npreg(bws=bw1,gradients=TRUE)
  
  gr <- gradients(kre1)
  
  dif1 <- outer(x[,1],x[,1],"-")
  dif2 <- outer(x[,2],x[,2],"-")
  dif <- sqrt(dif1^2 + dif2^2)
  
  kern <- function(x) (abs(x)<=1)*3/pi*(1-x^2)^2
  kern_deriv <- function(x) (abs(x)<=1)*12/pi*(1-x^2)
  
  lscv <- function(h) {
  der1 <- yt %*% t(rep(1, n)) *dif1/h^2*kern_deriv(dif/h)
  der2 <- yt %*% t(rep(1, n)) *dif2/h^2*kern_deriv(dif/h)
  den <- kern(dif/h)
  num <- rowSums(den) %*% t(rep(1, n)) - den
  gr0 <- cbind(colSums(der1 / num),colSums(der2 / num))
  return(sum(gr0^2)/n - 2*sum(gr0*gr)/n)
  }
  
  data_range <- max(apply(x,2,max)-apply(x,2,min))
  
  opt <- optimise(f = lscv, interval = c(data_range*0.02, data_range*0.5), tol = .Machine$double.eps)
  return(opt$minimum)
}

###############meanshift

##msOperator() implements one iteration of the regression mean shift algorithm at one location x_0
## The Epanechnikov kernel is used as the kernel g.
##Inputs: 
##x_0: current location; 
##x: covariate observations
##y: response observations
##z: kernel density estimates at x
##bw: bandwidth
msOperator <- function(x_0,x,y,z,bw){
  radi <- bw
  ne.id <- ball.id(x,x_0,radi)
  res <- c()
  if (length(ne.id)>0){
    ne <- x[ne.id,,drop=FALSE]
    n.ne <- length(ne.id)
    dif <- sweep(as.matrix(ne),2,as.matrix(x_0))
    dis <- sqrt(rowSums(dif^2))
    weig <- 1-(dis/bw)^2
    res <- colSums(weig*y[ne.id]/z[ne.id]*x[ne.id,,drop=FALSE])/sum(weig*y[ne.id]/z[ne.id])
  }else{
    res <- c(NA,NA)
  }
  res
}

##msAlgorithmUpdate() implements one iteration of the regression mean shift algorithm at multiple locations x0
##Inputs: 
##x0: current locations; 
##x: covariate observations
##y: response observations
##z: kernel density estimates at x
##bw: bandwidth
msAlgorithmUpdate <- function(x0,x,y,z,bw){
  M <- x0
  n <- nrow(x0)
  
  for(i in 1:n){
    M[i,] <- msOperator(x_0=x0[i,], x=x,y=y,z=z, bw=bw)
  }
  
  output <- M
  return( output )
}


#RegMS() implements the regression mean shift until convergence
##Inputs: 
##x: covariate observations
##yt: (transformed) response observations
##bw: bandwidth
##tol.stop: stopping threshold for the algorithm
##show.progress: indicator of whether or not to show the number of iterations
##show.graph: indicator of whether or not to plot the final results 
##background: indicator of whether or not to add a colored background in the plot
##cex.adj: indicator of whether or not to adjust the size of dots by the response variables in the plot
##plot.main: indicator of whether or not to add a title to the plot
##x1lim, x12lim: xy-range in the plot
RegMS <- function(x,yt,bw,tol.stop=1e-6,show.progress=FALSE,show.graph=FALSE,background=FALSE,cex.adj=FALSE,plot.main=TRUE,x1lim,x2lim){
  if (missing(x1lim)){
    x1lim=range(x[,1])
    x2lim=range(x[,2])
  } 
  
  if(plot.main==FALSE){
    main=""
  }else{
    main=paste("h=",round(bw,1),sep="")
  }
  
  if(background){
    x1grid=seq(x1lim[1],x1lim[2],length.out = 200)
    x2grid=seq(x2lim[1],x2lim[2],length.out = 200)
    
    fout <- outer(x1grid, x2grid, Vectorize(func))
    
    image(x1grid,x2grid,fout, col=terrain.colors(40), xlab="(a)", ylab="", main=paste("h=",round(bw,1),sep=""))
  }else if(show.graph){
    plot(x,xlim=x1lim,ylim=x2lim, type="n",main=main)
  }
  
  if(show.graph){
    if (cex.adj==TRUE){
      points(x,cex=yt)
    }else{
      points(x)
    }
  }
  
  dens <- kde(dat=x,h=bw,x=x)
  modes.old <- x
  
  stop.id <- rep(0,nrow(modes.old))
  k <- 1
  
  modes.final <- modes.old
  track.id <- 1:nrow(modes.old)
  while(sum(!stop.id)>0){
    if (show.progress){
      cat(k,",",sum(!stop.id),"\n")
    }
    modes.new <- msAlgorithmUpdate(modes.old,x,yt,dens,bw)
    if(show.graph){
      segments(modes.old[,1],modes.old[,2],modes.new[,1],modes.new[,2],col="red")
    }
    stop.id <- (sqrt(rowSums((modes.old - modes.new)^2))<tol.stop)
    modes.old <- modes.new[!stop.id,,drop=FALSE]
    if (sum(stop.id)>0){
      delete.id <- track.id[stop.id]
      track.id <- track.id[!stop.id]
      modes.final[delete.id,] <- modes.new[stop.id,]
    }
    k <- k+1
  }
  
  if (show.graph){
    points(modes.final[,1],modes.final[,2],col="blue",pch=4)
  }
  modes.final
}


##The following two functions distanceFunction() and connectedComponents() 
##are used to determine the connected components from the last points of all the mean shift trajactories.
##They are part of the archived R package MeanShift.
##For convenience, here they were directly copied from https://rdrr.io/cran/MeanShift/src/R/auxiliary.R

distanceFunction <-
  function (x, y) 
  {
    output <- sqrt(sum((x - y)^2))
    return(output)
  }

connectedComponents <- function( X, tol.epsilon=1e-3 ){
  
  N <- ncol( X )
  
  ## initialize components matrix
  C <- X
  
  ## initialize components vector
  labels <- vector( mode="integer", length=N )
  
  K <- 1 
  labels[1] <- 1
  C[,1] <- X[,1]
  
  # pb <- txtProgressBar( min=0, max=N, style=3 )
  
  ## efficient connected component algorithm
  for( n in 2:N ){
    
    assigned <- FALSE
    
    for( k in 1:K ){
      
      distance <- distanceFunction( X[,n], C[,k] )
      
      if( distance < tol.epsilon ){
        
        labels[n] <- k
        assigned <- TRUE
        break
        
      }
      
    }
    
    if( !assigned ){
      
      K <- K + 1
      labels[n] <- K
      C[,K] <- X[,n]
      
    }
    
    # setTxtProgressBar( pb, n )
    
  }
  
  C <- as.matrix( C[,1:K] )
  colnames( C ) <- paste( "mode", 1:K, sep="" )
  
  labels <- as.integer( labels )
  
  output <- list( components=C, labels=labels )
  
  # close( pb )
  
  message( "\nThe algorithm found ", as.character( K ),
           " clusters.\n")
  
  return( output )
  
}