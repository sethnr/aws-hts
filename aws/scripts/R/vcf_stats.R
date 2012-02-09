# use ggtools to read in VCF, calc fst across csome.
# having previously run RSetUp.R

library(GGtools)
gzpath <- "/Users/seth/Work/nextGen_snp_analysis/ngousso_emr/test_s_5_1.2sort.vcf.gz"
# chrom <- "2L"
# tabixcmd <- "~/Work/bin/tabix -c##"

vcftab <- TabixFile(gzpath, index = paste(gzpath, "tbi", sep="."))

gr2L <- GRanges(seqnames="0", IRanges(0,49364325))
gr2R <- GRanges(seqnames="1", IRanges(0,61545105))
gr3L <- GRanges(seqnames="2", IRanges(0,41963435))
gr3R <- GRanges(seqnames="3", IRanges(0,53200684))
grUNKN <- GRanges(seqnames="4", IRanges(0,42389979))
grX <- GRanges(seqnames="5", IRanges(0,24393108))
grY <- GRanges(seqnames="6", IRanges(0,237045))


grX.1<- GRanges(seqnames="5",IRanges(0,1e6))
#vcf2sm(gzpath, chrom, tabixcmd, nmetacol = 9, verbose = TRUE, gran=10000,metamax=100, makelocpref="chr")

smX <- vcf2sm(vcftab, gr=grX.1, nmetacol = 9L)
