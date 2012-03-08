# use ggtools to read in VCF, calc fst across csome.
# having previously run RSetUp.R

library(GGtools)
library(Rsamtools)
library(zoo)

wdpath <- "/Users/seth/Work/nextGen_snp_analysis/deepseq_fdr_3_9/";
setwd(wdpath)

vcfs <- getVcfs("FD0.*.vcf.gz$")
medDP <- getDPStats(vcfs)[[3]]
bafs <- getBafs(vcfs, minDP=(medDP/2), maxDP=(medDP * 2));


#plot rolling windows of 100snps, 5 snps step size
#plot(x=rollapply(bafs$POS,100,max, by=5), y=rollapply(abs(bafs$FD03_high.peaks - bafs$FD03_zero.peaks),100,mean, by=5), xlim=c(33002251,33661647))

vep <- getVep("./split/all.vep.sev.txt","./split/all.vep.txt")
data <- merge(vep, bafs, all.y=TRUE ,by="Uploaded_variation")
                                        #change factor to numeric values to numeric values
data$POS <- as.numeric(as.character(data$POS))

coding_data <- data[ which(!(data$Max_consequence %in% c("UPSTREAM","DOWNSTREAM","INTRONIC"))),]
#remove any silent mutations:
nonsyn_data <- data[ which((data$Max_consequence %in% c("NON_SYNONYMOUS_CODING","FRAMESHIFT_CODING","SPLICE_SITE","ESSENTIAL_SPLICE_SITE","STOP_GAINED","STOP_LOST"))),]

# bind contrast to DF
#coding_data <- cbind(coding_data,Contrast=abs(coding_data$FD03_high.peaks - coding_data$FD03_zero.peaks))
qplot(POS,y=Contrast,data=coding_data,xlim=c(33002251,33661647),color=Max_consequence)


genesOnly <- getGenesBM()

#nb: requires full getGenesBM
go <- unique(as.data.frame(cbind(genes$ensembl_transcript_id,genes$go_accession,genes$name_1006)))
ipr <- unique(as.data.frame(cbind(genes$ensembl_transcript_id,genes$interpro,genes$interpro_description)))




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
