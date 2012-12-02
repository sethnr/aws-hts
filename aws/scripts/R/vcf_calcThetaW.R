source("http://bioconductor.org/biocLite.R")
biocLite("ShortRead")
biocLite("snpStats");
biocLite("GGtools");

library(ShortRead)
library(snpStats)
library(GGtools)

source("~/Work/EC2/aws-hts/aws/scripts/R/lib/VCFs.R");
filename <- "~/Scratch/mergeTest/combined-vars-only.vcf.gz";

vcf <- getVCF(filename)

vcf <- getBafsMultiFile(vcf)

calcThetaVcf(vcf=vcf)