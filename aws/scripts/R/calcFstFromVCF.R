library(GGtools);
library(SnpStats);
library(ShortRead);
getVcfParam <- function(chr, start, end) {
  param <- ScanVcfParam(which = GRanges(chr, IRanges(start, end)))
  param
}
#params = getVcfParam("4",6780000, 6790000);
#gr <- GRanges(chr, IRanges(start, end))
#vcfSM <- getSM(filename, grange=gr);

source("~/Work/EC2/aws-hts/aws/scripts/R/lib/VCFs.R");

filename <- "~/Scratch/mergeTest/combined.vars.vcf.gz";
vcf <- getVCF(filename);
vcf <- getBafsMultiFile(vcf)
bafs <- data.frame(vcf[[1]]$CHR, vcf[[1]]$POS, vcf[[1]]$GENO$BAF)

start <- 6780000;
end <- 6790000;
chr <- 4
#calcThetaVcf(vcf = vcf, noIndivs = 10);
calcThetaVcf(vcf = vcf, chr=chr, start=start, end=end, noIndivs = 10);

