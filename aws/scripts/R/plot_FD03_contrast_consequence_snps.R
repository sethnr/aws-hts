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

wdpath <- "/Users/seth/Work/nextGen_snp_analysis/deepseq_fdr_3_9/manu_peaks";
#wdpath <- "/Users/seth/Work/nextGen_snp_analysis/deepseq_fdr_3_9/APL1s";
setwd(wdpath)

vcfs <- getVcfs("FD0.*3LpeakF.*.vcf.gz$")
# get median depth
medDP <- getDPStats(vcfs)[[3]]
bafs <- getBafs(vcfs, minDP=(medDP/2), maxDP=(medDP * 2));

# add [high - low] as contrast value
bafs <- cbind(bafs,Contrast03=abs(bafs$FD03_3LpeakF_high - bafs$FD03_3LpeakF_zero))
bafs <- cbind(bafs,Contrast09=abs(bafs$FD09_3LpeakF_high - bafs$FD09_3LpeakF_zero))


#plot rolling windows of 100snps, 5 snps step size
#plot(x=rollapply(bafs$POS,100,max, by=5), y=rollapply(abs(bafs$FD03_high.peaks - bafs$FD03_zero.peaks),100,mean, by=5), xlim=c(33002251,33661647))

vep <- getVep("./FDs_3LpeakF_all.vep","./FDs_3LpeakF_all.sev.vep")
data <- merge(vep, bafs, all.y=TRUE ,by="Uploaded_variation")
                                        #change factor to umeric values to numeric values

oldNames <- names(data)
oldNames[3] <- "POS"
names(data) <- oldNames
data$POS <- as.numeric(as.character(data$POS))


# add col if trans_consequence is max_consequence
isMainCons <- apply(cbind(as.character(data$Max_consequence),as.character(data$Trans_consequence)), 1,FUN=function(x) {grepl(x[1],x[2])})
data <- cbind(data,isMainCons)
data <- data[ which(data$isMainCons ==1),]
data <- cbind(data,Contrast03=abs(coding_data$FD03_high.peaks - coding_data$FD03_zero.peaks))
data <- cbind(data,Contrast09=abs(coding_data$FD09_high.peaks - coding_data$FD09_zero.peaks))



coding_data <- data[ which(!(data$Max_consequence %in% c("UPSTREAM","DOWNSTREAM","INTRONIC","INTERGENIC"))),]
#remove any silent mutations:
nonsyn_data <- data[ which((data$Max_consequence %in% c("NON_SYNONYMOUS_CODING","FRAMESHIFT_CODING","SPLICE_SITE","ESSENTIAL_SPLICE_SITE","STOP_GAINED","STOP_LOST"))),]

# bind contrast to DF
#coding_data <- cbind(coding_data,Contrast=abs(coding_data$FD03_high.peaks - coding_data$FD03_zero.peaks))
#qplot(POS,y=Contrast,data=coding_data,xlim=c(33002251,33661647),color=Max_consequence)
qplot(POS,y=Contrast03,data=data,color=Max_consequence,xlim=c(33002251,33661647))
dev.print(png,file="FD03_3Lpeak_contrast.png",width=1200,height=800);
qplot(POS,y=Contrast09,data=data,color=Max_consequence,xlim=c(33002251,33661647))
dev.print(png,file="FD09_3Lpeak_contrast.png",width=1200,height=800);

qplot(POS,y=Contrast03,data=coding_data,color=Max_consequence,xlim=c(33002251,33661647))
dev.print(png,file="FD03_3Lpeak_coding_contrast.png",width=1200,height=800);
qplot(POS,y=Contrast09,data=coding_data,color=Max_consequence,xlim=c(33002251,33661647))
dev.print(png,file="FD09_3Lpeak_coding_contrast.png",width=1200,height=800);


