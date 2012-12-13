library(reshape2)
library(ggplot2)

plotThetaChr <- function(thetas, chr, filename) {
  thetasChr<-thetas[which(thetas$chr==chr),]
  thetas_m <- melt(thetasChr, id=c("chr","start","end"))
  if(exists("file")) { png(paste(filename,chr,"png",sep="."))}
  qplot(start,value,main=paste(chr," theta cf"),data=thetas_m,color=variable)
  if(exists("file")) { dev.off() }
  }

setwd("~/Scratch/relTheta/Fd09/")
infile <- "Fd09.theta"
thetas <- read.table(infile,header=T)

relThetas <- thetas[,4:6] / thetas[,7]
colnames(relThetas) <- paste(colnames(relThetas),"_r", sep="")
relThetas <- cbind(thetas[,1:3],relThetas)
write.table(relThetas,file="Fd09.theta.rel",col.names=T, row.names=F, quote=F);

for (chr in (unique(thetas$CHR))) {
  plotThetaChr(relThetas,chr,filename);
}

thetas.scale <- t(scale(t(thetas[,4:6]),scale=thetas[,7]))
colnames(thetas.scale) <- paste(colnames(thetas.scale),"_s", sep="")
thetas.scale<-cbind(thetas[,1:3],thetas.scale)
write.table(thetas.scale,file="Fd09.theta.sca",col.names=T, row.names=F, quote=F);
plotChr(thetas.scale,"2L");

