library(reshape2)
library(ggplot2)

plotThetaChr <- function(thetas, chr, file) {
  thetasChr<-thetas[which(thetas$chr==chr),]
  thetas_m <- melt(thetasChr, id=c("chr","start","end"))
  png(paste(filename,chr,"png",sep="."))
  qplot(start,value,main=paste(chr," theta cf"),data=thetas_m,color=variable)
  dev.off()
  }

setwd("~/Scratch/relTheta_hp/")


runThetaGraphs <- function(infile) {
  thetas <- read.table(infile,header=T)
  
  relThetas <- thetas[,4:6] / thetas[,7]
  colnames(relThetas) <- paste(colnames(relThetas),"_r", sep="")
  relThetas <- cbind(thetas[,1:3],relThetas)
  relFile <- paste(infile,".rel",sep="")
  write.table(relThetas,file=relFile,col.names=T, row.names=F, quote=F);
  
  for (chr in (unique(thetas$CHR))) {
    plotThetaChr(relThetas,chr,filename);
  }
  
  thetas.scale <- t(scale(t(thetas[,4:6]),scale=thetas[,7]))
  colnames(thetas.scale) <- paste(colnames(thetas.scale),"_s", sep="")
  thetas.scale<-cbind(thetas[,1:3],thetas.scale)
  scaFile <- paste(infile,".sca",sep="")
  write.table(thetas.scale,file=scaFile,col.names=T, row.names=F, quote=F);
 # plotThetasChr(thetas.scale,"2L",scaFile);
  
  
  #find Sds for each 3 pools
  theta_sd <- data.frame(chr=thetas$chr, start=thetas$start, end=thetas$end,theta.sd=apply(thetas.scale[,4:6],1,sd))
  # make granges, segment and change to data.frame
  theta_sd_gr <- GRanges(seqnames=Rle(theta_sd$chr), ranges=IRanges(theta_sd$start, end=theta_sd$end), depth=theta_sd$theta.sd)
  theta_sd_segs <- fastseg(theta_sd_gr)
  theta_sd_segs <- as.data.frame(theta_sd_segs)[,c("seqnames","start","end","seg.mean")]
  
  #get 3rd q for topSegs
  q3 <- summary(theta_sd$theta.sd)[[5]]
  
  # change names and write tables
  #sd
  sdFile <- paste(infile,".sd",sep="")
  colnames(theta_sd) <- c("chr","start","end",paste(infile,".sd",sep=""))
  write.table(theta_sd,file=sdFile,col.names=T,row.names=F,quote=F)
  #segs
  colnames(theta_sd_segs) <- c("chr","start","end",paste(infile,".sd_seg",sep=""))
  sdSegFile <- paste(infile,".sd.segs",sep="")
  write.table(theta_sd_segs,file=sdSegFile,col.names=T,row.names=F,quote=F)
  #topsegs
  write.table(theta_sd_segs[which(theta_sd_segs[4] >= (q3 *2)),],file=paste(infile,".sd.Q3_2",sep=""),col.names=T,row.names=F,quote=F)
  write.table(theta_sd_segs[which(theta_sd_segs[4] >= (q3 *3)),],file=paste(infile,".sd.Q3_3",sep=""),col.names=T,row.names=F,quote=F)
  write.table(theta_sd_segs[which(theta_sd_segs[4] >= (q3 *4)),],file=paste(infile,".sd.Q3_4",sep=""),col.names=T,row.names=F,quote=F)
}


runThetaGraphs(infile="Fd03.theta")
runThetaGraphs(infile="Fd09.theta")

