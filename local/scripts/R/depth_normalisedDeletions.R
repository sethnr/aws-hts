#get relative depth drops:
setwd("/Users/seth/Scratch/dmfa_imls")
depthG.f <-"Pfin5G/Pfin5G.vcfD"
depthP.f <-"Pfin5P/Pfin5P.vcfD"

depthP <- read.table(depthP.f, header=T)
depthG <- read.table(depthG.f, header=T)



depth <- merge(depthG[,1:4],depthP[,1:4])
quants <- rowSums(depth[,3:6])/4
quants <- quants[order(quants)]

#do QQ plots
if(1==2) {
  sample <- sample(dim(depth)[[1]],10000)
  plot(qqplot(depth[,3],depth[,4]))
  plot(qqplot(depth[,5],depth[,6]))
  plot(qqplot(rowSums(depth[sample,3:4]),rowSums(depth[sample,5:6])))
}

for (i in 3:6) {
  order <- order(depth[,i])
  depth[order,i] <- quants
}

depth <- depth[order(depth[,1],depth[,2]),]
depth <- cbind(depth, Pfin5G.total=rowSums(depth[3:4]),Pfin5P.total=rowSums(depth[5:6]))
#depth.scale.file <- paste(depth.file,".s", sep="")
depth <- cbind(depth,relDep.high.G =depth[,3]/(depth[,7]/2),relDep.high.P=depth[,5]/(depth[,8]/2))

depth.2LA <- subset(depth,chr=="2L" & pos>=21199039 & pos <=42127751)
sample <- sample(1000000)
plot(density(depth[sample,3]),xlim=c(0,150),main="depth density: 1m sample genome, 100k sample 2La" )
lines(density(depth[sample,4]))
lines(density(depth[sample,5]))
lines(density(depth[sample,6]))
sample <- sample(100000)
lines(density(depth.2LA[sample,3]),col="red")
lines(density(depth.2LA[sample,4]),col="red")
lines(density(depth.2LA[sample,5]),col="red")
lines(density(depth.2LA[sample,6]),col="red")
                  
#depth <- head(depth,1000)
#shift each one by mean distance from half total (i.e. normalize overall depth)
#depth.scale <- scale(depth[3:4], center=(colMeans(depth[,3:4])-mean(depth[,5]/2)),scale=F)
write.table(depth, file="Pfin5.depth.qnorm", col.names=T, row.names=F, quote=F)
