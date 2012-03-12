# use ggtools to read in VCF, calc fst across csome.
# having previously run RSetUp.R

library(GGtools)
library(Rsamtools)
library(zoo)
library(ggplot2)


setwd("~/Work/EC2/aws-hts/aws/");
source("scripts/R/lib/VCFs.R");
source("scripts/R/lib/Annotation.R");
source("scripts/R/lib/CNVs.R");

wdpath <- "/Users/seth/Work/nextGen_snp_analysis/deepseq_fdr_3_9/";
setwd(wdpath)

vcfs <- getVcfs("FD0.*.vcf.gz$")
# get median depth
medDP <- getDPStats(vcfs)[[3]]
bafs <- getBafs(vcfs, minDP=(medDP/2), maxDP=(medDP * 2));

# add [high - low] as contrast value
bafs <- cbind(bafs,Contrast=abs(bafs$FD03_high.peaks - bafs$FD03_zero.peaks))


#plot rolling windows of 100snps, 5 snps step size
#plot(x=rollapply(bafs$POS,100,max, by=5), y=rollapply(abs(bafs$FD03_high.peaks - bafs$FD03_zero.peaks),100,mean, by=5), xlim=c(33002251,33661647))

vep <- getVep("./split/all.vep.txt","./split/all.vep.sev.txt")
data <- merge(vep, bafs, all.y=TRUE ,by="Uploaded_variation")
                                        #change factor to numeric values to numeric values
data$POS <- as.numeric(as.character(data$POS))

# add col if trans_consequence is max_consequence
isMainCons <- apply(cbind(as.character(data$Max_consequence),as.character(data$Trans_consequence)), 1,FUN=function(x) {grepl(x[1],x[2])})
data <- cbind(data,isMainCons)
coding_data <- data[ which(!(data$Max_consequence %in% c("UPSTREAM","DOWNSTREAM","INTRONIC","INTERGENIC"))),]
#remove any silent mutations:
nonsyn_data <- data[ which((data$Max_consequence %in% c("NON_SYNONYMOUS_CODING","FRAMESHIFT_CODING","SPLICE_SITE","ESSENTIAL_SPLICE_SITE","STOP_GAINED","STOP_LOST"))),]

# bind contrast to DF
#coding_data <- cbind(coding_data,Contrast=abs(coding_data$FD03_high.peaks - coding_data$FD03_zero.peaks))
qplot(POS,y=Contrast,data=coding_data,xlim=c(33002251,33661647),color=Max_consequence)


genes <- getGenesBM()
genesGO <- getGenesGO()
genesIpr <- getGenesIpr()

# get non-synonymous vars & merge with gene GO terms
diffNonsyn <- merge(nonsyn_data[which(nonsyn_data$Contrast >= 0.6),c("Uploaded_variation","CHR","POS","Gene","Feature","Trans_consequence","Max_consequence","Contrast")], genes, all.x=TRUE, by.x="Feature", by.y="ensembl_transcript_id")




##########
##
##   DN:DS stuff
##
##########


dnds <- getDNDS(vep)
q3 <- summary(dnds$dnds)[[5]]
maxDnds <- dnds[which(dnds$dnds >= q3),]
maxDnds <- merge(maxDnds, genesOnly[which(genesOnly$ensembl_transcript_id %in% maxDnds$Feature),c("ensembl_transcript_id","description")], by.x="Feature",by.y="ensembl_transcript_id")

q1 <- summary(dnds$dnds)[[2]]
minDnds <- dnds[which(dnds$dnds <= q1),]
minDnds <- merge(minDnds, genesOnly[which(genesOnly$ensembl_transcript_id %in% minDnds$Feature),c("ensembl_transcript_id","description")], by.x="Feature",by.y="ensembl_transcript_id")
rm(q1,q3)


dps <- getDP(vcfs)
dpsM <- melt(dps, id=c("CHR","POS","Uploaded_variation"))
qplot(POS,y=value,data=dpsM,color=variable, xlim=c(33002251,33661647))


