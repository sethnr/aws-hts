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

#nas <- is.na(vcfSM[1:3,])
#nas <- apply(is.na(vcfSM[1:3,]),2,sum)
#vcfSM[1:3,which(nas <= 3)]

noNa03 <- which(apply(is.na(vcfSM[1:3,]),2,sum) < 3)
noNa09 <- which(apply(is.na(vcfSM[4:6,]),2,sum) < 3)
ld.vcf.03 <- ld(vcfSM[1:3,noNa03],depth=5000, stats=c("R.squared"))
ld.vcf.09 <- ld(vcfSM[4:6,noNa09],depth=5000, stats=c("R.squared"))

tags03 <- colSums(ld.vcf.03, na.rm = TRUE, dims = 1)
tags09 <- colSums(ld.vcf.09, na.rm = TRUE, dims = 1)

names <- vcfSM@dimnames[[2]]
posns <- as.numeric(sub("chr3L:", "", names))

plot(posns[noNa03], tags03)
dev.print(png,file="3LPeakFlank_R.squared_noNA_FD03.png",width=1200,height=800)

plot(posns[noNa09], tags09)
dev.print(png,file="3LPeakFlank_R.squared_noNA_FD09.png",width=1200,height=800)



ld.vcf <- ld(vcfSM,depth=5000, stats=c("D.prime","R.squared"))
tagsAll <- colSums(ld.vcf$R.squared, na.rm = TRUE, dims = 1)

window <- 10000
step <- 2000
maxes <- rollapply(tagsAll,window,max,by=step)
indices <- which(tagsAll %in% maxes)
plot(indices, unique(maxes))
