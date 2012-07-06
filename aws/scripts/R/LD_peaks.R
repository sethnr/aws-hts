# use ggtools to read in VCF, calc fst across csome.
# having previously run RSetUp.R

library(GGtools)
library(Rsamtools)
library(hexbin)


setwd("~/Work/EC2/aws-hts/aws/");
source("scripts/R/lib/VCFs.R");
source("scripts/R/lib/Annotation.R");
source("scripts/R/lib/CNVs.R");
source("scripts/R/getAgamGRanges.R");
wdpath <- "/Users/seth/Work/nextGen_snp_analysis/deepseq_fdr_3_9/manu_peaks";
setwd(wdpath)


pattern <- "FDs.*.vcf.gz$"
vcfFiles <- list.files(path = ".", pattern = pattern)
#  vcfFiles
#  [1] "FDs_3LpeakF_all.phase.vcf.gz"

gr3L <- GRanges(seqnames="3L", IRanges(0,41963435))

spectrum <- rainbow(10, start=0, end=1/6)[10:1]
filename <- vcfFiles[[1]]
vcftab <- TabixFile(filename, index = paste(filename, "tbi", sep="."));
vcfSM <- vcf2sm(vcftab, gr=gr3L, nmetacol=9L)


ld.vcf <- ld(vcfSM,depth=1000, stats=c("D.prime","R.squared"))
image(ld.vcf$D.prime, lwd=0, cuts=9, col.regions=spectrum, colorkey=TRUE)
dev.print(png,file="3LPeakFlank_D.prime_1k.png",width=1200,height=800)

ld.vcf <- ld(vcfSM,depth=2000, stats=c("D.prime","R.squared"))
image(ld.vcf$D.prime, lwd=0, cuts=9, col.regions=spectrum, colorkey=TRUE)
dev.print(png,file="3LPeakFlank_D.prime_10k.png",width=1200,height=800)

m <- ld.vcf
ld.vcf.xy <- data.frame(r=as.vector(row(m)), c=as.vector(col(m)), v=as.vector(m))
rm(m)

snpNames <- ld.vcf$D.prime@Dimnames[[1]]
#snpInfo <- as.data.frame(strsplit(snpNames,"[:]"))

os <- vector();
chr <- vector()
 for(i in 1:length(snpNames)) {
     vls <- strsplit(as.character(snpNames[[i]]),"[:]")[[1]]
     pos <- c(pos,as.numeric(vls[2]))
     chr <- c(chr,vls[1])
   }

# get posn of SNP of interest in vector
posn <- which(pos==33073484,arr.ind=TRUE)

w <- 100
diags <- vector("list", w)
for (i in 1:w) diags[[i]] <- pos[(i+1):posn] - pos[1:(posn-i)]
dist <- bandSparse(posn, k=1:w, diagonals=diags)

distance <- dist@x
D.prime <- ld.vcf$D.prime@x
R.squared <- ld.vcf$R.squared@x


# look at decay of R^2 around a point
lr100 <- c((posn-w):(posn-1), (posn+1):(posn+w))
D.prime <- ld(vcfSM[,posn], vcfSM[,lr100], stats="D.prime")
where <- pos[lr100]

D.prime[is.na(D.prime)] <- 0 
plot(where, D.prime)
lines(where, smooth(D.prime))


r2 <- forceSymmetric(ld.vcf$R.squared[lr100, lr100])
D <- as.dist(1-r2)
hc <- hclust(D, method="complete")
par(cex=0.5)
plot(hc)




ld.ceph <- ld(ceph.1mb, stats=c("D.prime", "R.squared"), depth=100)
diags <- vector("list", 100)
for (i in 1:100) diags[[i]] <- pos[(i+1):603] - pos[1:(603-i)]
dist <- bandSparse(603, k=1:100, diagonals=diags)

distance <- dist@x
D.prime <- ld.ceph$D.prime@x
R.squared <- ld.ceph$R.squared@x





# make block-means from array
# x <- matrix(1:100, 10)  # create data
x <- ld.vcf$D.prime
b<-100
md <- dim(x)[1]
newD <- round(md/b)
rmean <- matrix(0,newD,newD)  # result matrix
for (i in 1:newD){
     for (j in i:newD){
         bj <-  c(-b+1,0) + b * j
         bi <- c(-b+1,0) + b * i
         if(bj[[2]] > md) bj[[2]]<-md
         if(bi[[2]] > md) bi[[2]]<-md
         block <- x[bi[[1]]:bi[[2]],bj[[1]]:bj[[2]]]
         block[isNA(block)] <- 0
         rmean[i, j] <- sum(block)
       }
     write(paste("block ",i," ",j," -> ",bi," ",bj),stderr())
        
 }
x
rmean

image(rmean,breaks=seq(0,10000,length=(l+1)),col=rgb(red=1,green=seq(1,0,length=l),blue=seq(1,0,length=l)),axes=FALSE);
axis(1, at=seq(0,1,length=6), labels = seq(32002253, 34661703, length=6))


ld.vcf <- ld(vcfSM,depth=5000, stats=c("D.prime","R.squared"))
tagsX <- colSums(ld.vcf$R.squared, na.rm = TRUE, dims = 1)
plot(tagsX)

meandist <- ((34661703 - 32002253) / (90765/5000))/2
breaks <- seq((32000000 + meandist),(35000000 -meandist),length=20)
hcs <- list()

for (i in 2:length(breaks)) {
  #plot cluster to find tag SNPs, 

  write(paste("getting range: ",breaks[(i-1)]-meandist," - ",breaks[i]+meandist),stderr())      
  grange <- GRanges(seqnames="3L", IRanges(breaks[(i-1)],breaks[i]))
  vcfSM <- vcf2sm(vcftab, gr=grange, nmetacol=9L)
  ld.vcf <- ld(vcfSM,depth=5000, stats=c("D.prime","R.squared"))
  r2<-forceSymmetric(ld.vcf$R.squared)
  r2[is.na(r2)] <- 0
  D <- as.dist(1-r2)
  hc <- hclust(D, method="ave")
  par(cex=0.5)
  hcs <- c(hcs, hc)
  plot(hc)
}


# try rollapply (which.max(x)) to find local peaks, get array, then uniq() it. See how long list is & positioning.
vcfSM <- vcf2sm(vcftab, gr=gr3L, nmetacol=9L)
ld.vcf <- ld(vcfSM,depth=5000, stats=c("D.prime","R.squared"))

tagsX <- colSums(ld.vcf$R.squared, na.rm = TRUE, dims = 1)


allTags <- vector()
limit <- 0.9
for (i in 1:length(hcs)) {
  allTags <- unique(c(allTags, names(which(cutree(hc, h=limit)==1))))
}
indices <- which(ld.vcf$R.squared@Dimnames[[1]] %in% allTags, arr.ind=TRUE)
plot(tagsX[indices])


window <- 10000
step <- 2000
maxes <- rollapply(tagsX,window,max,by=step)
indices <- which(tagsX %in% maxes)
plot(indices, unique(maxes))
