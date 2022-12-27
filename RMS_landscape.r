##This R code reproduces Figure 4.3

source("RMS_funs.r")
#############
library(np)
library(rgl)
library(colorRamps)

#status <- "WT"
status <- "Q61L"

data.sample <- read.table(paste("proteins/", status,".txt",sep=""),header =TRUE)
#data available at https://dx.doi.org/10.21227/331n-7019.
n <- nrow(data.sample)

x <- data.sample[,1:2]
y <- data.sample[,3]
yt <- abs(y)

bw <- gbbw(x,yt)
bw 

modes.final <- RegMS(as.matrix(x),yt,bw)
output <- connectedComponents( X=t(modes.final))

symb <- c(1,6,8)

layout(matrix(1:2,ncol=2), width = c(10,1),height = c(1,1))
par(mar=c(4,4,1,0), oma = c(1, 1, 1, 1))

red.threshold <- -6000
blue.threshold <- -6600
id.middle <- (1:n)%in%which(y < red.threshold & y > blue.threshold)
id.blue <- (1:n)%in%which(y <= blue.threshold)
id.red <- (1:n)%in%which(y >= red.threshold)

palette <- rev(matlab.like2(1000))
Col <- palette[as.numeric(cut(yt[id.middle],breaks = 1000))]
plot(x, type="n", col=Col, cex=0.8,bty="l",xlab="PC1", ylab="PC2",main=status)
points(x[id.red, ],pch=symb[output$labels[id.red]],col="#BF0000",cex=0.8)
points(x[id.middle, ], pch=symb[output$labels[id.middle]], col=Col, cex=0.8)
points(x[id.blue, ],pch=symb[output$labels[id.blue]],col="#0000BF",cex=0.8)

points(modes.final[,1],modes.final[,2],pch=4,cex=2)

par(mar=c(0,0,0,0))
legend_image <- as.raster(matrix(palette, ncol=1))
plot(c(0,0.4),c(0,1),type = 'n', axes = F,xlab = '', ylab = '')
text(x=0.27, y = seq(0.10,1,l=5), labels = round(seq(blue.threshold,red.threshold,l=5)))
rasterImage(legend_image, 0, 0.10, 0.1,1)






