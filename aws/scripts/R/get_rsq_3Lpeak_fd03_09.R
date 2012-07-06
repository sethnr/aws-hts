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

filename <- vcfFiles[[1]]
vcftab <- TabixFile(filename, index = paste(filename, "tbi", sep="."));
vcfSM <- vcf2sm(vcftab, gr=gr3L, nmetacol=9L)

ld.03 <- ld(vcfSM[1:3,],depth=5000, stats=c("R.squared"))
ld.09 <- ld(vcfSM[4:6,],depth=5000, stats=c("R.squared"))

rsq03 <- colSums(ld.03, na.rm = TRUE, dims = 1)
rsq09 <- colSums(ld.09, na.rm = TRUE, dims = 1)
posns <- sapply(strsplit(dimnames(vcfSM)[[2]],":",fixed=TRUE),"[[",2)

plot(posns,rsq03)
dev.print(png,file="3LPeakFlank_R.squared_FD03.png",width=1200,height=800)
plot(posns,rsq09)
dev.print(png,file="3LPeakFlank_R.squared_FD09.png",width=1200,height=800)


window <- 10000
step <- 2000

maxes <- rollapply(rsq03,window,max,by=step)
plot(posns[which(rsq03 %in% maxes)], unique(maxes))
dev.print(png,file="3LPeakFlank_R.sq_tags_FD03.png",width=1200,height=800)

maxes <- rollapply(rsq09,window,max,by=step)
plot(posns[which(rsq09 %in% maxes)], unique(maxes))
dev.print(png,file="3LPeakFlank_R.sq_tags_FD09.png",width=1200,height=800)
