##This R code reproduces Figure 4.4

source("RMS_funs.r")
#############
library(SPODT)

data(dataMALARIA)
x <- dataMALARIA[,2:3]
names(x) <- c("Longitude", "Latitude")
y <- dataMALARIA[,4]
yt <- y + 0.1

bw <- gbbw(x,yt)
bw

par(mfrow=c(1,2))
par(mar=c(5,5,1,0), oma = c(1, 1, 1, 1))
modes.final <- RegMS(as.matrix(x),yt,bw)
output <- connectedComponents( X=t(modes.final))
colors <- c("red", "black", "green")
plot(x,col=colors[output$labels], cex=yt+0.3,pch=19)


##############

library(tree)
fit = tree(y~x[,1]+x[,2])

plot(x,cex=yt+0.3,pch=19)
partition.tree(fit,add=TRUE)
