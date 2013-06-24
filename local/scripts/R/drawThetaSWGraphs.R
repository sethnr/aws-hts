library(reshape2)
library(ggplot2)
library(GenomicRanges)
library(fastseg)

plotThetaChr <- function(thetas, chr, file) {
  thetasChr<-thetas[which(thetas$chr==chr),]
  thetas_m <- melt(thetasChr, id=c("chr","start","end"))
  png(paste(filename,chr,"png",sep="."))
  qplot(start,value,main=paste(chr," theta cf"),data=thetas_m,color=variable)
  dev.off()
  }

setwd("~/Scratch/relTheta_hp/")


runThetaGraphs <- function(infile) {
  thetas <- read.table(infile,header=TRUE)
  if(!colnames(thetas)[[1]] == "chr") {
    thetas <- read.table(infile,header=FALSE)
    colnames(thetas) <- c("chr","start","end",paste(infile, c("high","low","zero","total"), sep="."))
    }
  relThetas <- thetas[,4:6] / thetas[,7]
  colnames(relThetas) <- paste(colnames(relThetas),"_r", sep="")
  relThetas <- cbind(thetas[,1:3],relThetas)
  relFile <- paste(infile,".rel",sep="")
  write.table(relThetas,file=relFile,col.names=TRUE, row.names=FALSE, quote=FALSE);
  
#  for (chr in (unique(thetas$CHR))) {
#    plotThetaChr(relThetas,chr,filename);
#  }
  
  thetas.scale <- t(scale(t(thetas[,4:6]),scale=thetas[,7]))
  colnames(thetas.scale) <- paste(colnames(thetas.scale),"_s", sep="")
  thetas.scale<-cbind(thetas[,1:3],thetas.scale)
  scaFile <- paste(infile,".sca",sep="")
  write.table(thetas.scale,file=scaFile, col.names=TRUE, row.names=FALSE, quote=FALSE);
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
  write.table(theta_sd,file=sdFile,col.names=TRUE,row.names=FALSE,quote=FALSE)
  #segs
  colnames(theta_sd_segs) <- c("chr","start","end",paste(infile,".sd_seg",sep=""))
  sdSegFile <- paste(infile,".sd.segs",sep="")
  write.table(theta_sd_segs,file=sdSegFile,col.names=TRUE,row.names=FALSE,quote=FALSE)
  #topsegs
  write.table(theta_sd_segs[which(theta_sd_segs[4] >= (q3 *2)),],file=paste(infile,".sd.Q3_2",sep=""),col.names=TRUE,row.names=FALSE,quote=FALSE)
  write.table(theta_sd_segs[which(theta_sd_segs[4] >= (q3 *3)),],file=paste(infile,".sd.Q3_3",sep=""),col.names=TRUE,row.names=FALSE,quote=FALSE)
  write.table(theta_sd_segs[which(theta_sd_segs[4] >= (q3 *4)),],file=paste(infile,".sd.Q3_4",sep=""),col.names=TRUE,row.names=FALSE,quote=FALSE)
  theta_sd
}


theta_sd <- runThetaGraphs(infile="Fd03_a.theta")
runThetaGraphs(infile="Fd09_a.theta")





##uniquify perms##
thetaP <- read.table("~/Scratch/FD03/thetaP_a/Fd03_P_a.theta")

theta_100 <- data.frame(row.names=rownames(theta))
for (start in (unique(thetaP$start))) {
  posSample <- thetaP[which(thetaP$start == start),]
  noSamples <- dim(posSample)[[1]]
  write(paste(start, noSamples),stderr())
  posSample <- posSample[sample(1:noSamples,100),]
  thetaP_100 <- rbind(thetaP_100,posSample)
}
thetaP_100.scale <- t(scale(t(thetaP_100[,4:6]),scale=thetaP_100[,7]))
colnames(thetaP_100.scale) <- paste(colnames(thetaP_100.scale),"_s", sep="")
thetaP_100.scale<-cbind(thetaP_100[,1:3],thetaP_100.scale)

thetaP_100_sd <- data.frame(chr=thetaP_100$chr, 
                       start=thetaP_100$start, 
                       end=thetaP_100$end,
                       thetaP.sd=apply(thetaP_100.scale[,4:6],1,sd))

qplot(thetaP_100_sd$start+5000,thetaP_100_sd[,4])
plot(density(thetaP_100_sd$theta.sd))
