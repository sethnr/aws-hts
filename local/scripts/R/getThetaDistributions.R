library(reshape2)
library(ggplot2)
library(GenomicRanges)
library(fastseg)


getThetaSd <- function(thetas) {
  #change thetas to numeric and transpose
  rawThetas <- rbind(as.numeric(thetas[,4]),as.numeric(thetas[,5]),as.numeric(thetas[,6]))
  thetas.scale <- t(scale(rawThetas,scale=as.numeric(thetas[,7])))
#  thetas.scale<-cbind(thetas[,1:3],thetas.scale)
#  theta_sd <- data.frame(chr=thetas$chr, start=thetas$start, end=thetas$end,theta.sd=apply(thetas.scale[,4:6],1,sd))
  theta.sd=apply(thetas.scale,1,sd)  
  theta.sd
}


#fin <- file("stdin")

#fin <- file("~/Scratch/FD03/thetaP_a/Fd03_P_a.theta")
#open(fin)
#fout <- file("~/Scratch/FD03/thetaP_a/Fd03_P_a.theta.sd", "w")
#flocout <- file("~/Scratch/FD03/thetaP_a/Fd03_P_a.theta.sd.loc", "w")

fin <- file("~/Scratch/FD03/thetaP_a/Fd03_P_a.theta.rem")
open(fin)
fout <- file("~/Scratch/FD03/thetaP_a/Fd03_P_a.theta.rem.sd", "w")
flocout <- file("~/Scratch/FD03/thetaP_a/Fd03_P_a.theta.rem.sd.loc", "w")
sd <- vector()
rm(locMatrix)
while(length(line <- readLines(fin,n=1)) > 0) {
  lineSplit <- strsplit(line, "\t")[[1]]
  if(exists("locMatrix")) {
    if(locMatrix[1,2] == lineSplit[[2]]) {
      locMatrix <- rbind(locMatrix, lineSplit)
    }
    if(locMatrix[1,2] != lineSplit[[2]]) {
      locSd <- getThetaSd(locMatrix)
#      sd <- append(sd,locSd)
      write(locSd, fout, ncolumns=1)
      outline <- c(locMatrix[1,1:3],quantile(locSd,c(0.50,0.95,0.99)))
      write(paste(outline,sep="\t"), flocout, ncolumns=length(outline))
      locMatrix <- matrix(lineSplit, ncol=7)
    }
  }
  if(!exists("locMatrix")) {
    locMatrix <- matrix(lineSplit, ncol=7)    
  }
}
locSd <- getThetaSd(locMatrix)
#sd <- append(sd,locSd)
write(fout, paste(locSd,"\n"))
outline <- c(locMatrix[1,1:3],quantile(locSd,c(0.50,0.95,0.99)))
write(paste(outline,sep="\t"), flocout)
close(fout)
close(flocout)
close(fin)


sds <- scan(file="/Users/seth/Scratch/FD03/thetaP_a/Fd03_P_a.theta.c.sd", what=numeric())
plot(density(sds))
quantile(sds,probs=c(0.5, 0.95, 0.99, 0.999))