library(zoo)
library(Rsamtools)
library(ggplot2)

setwd("~/Work/EC2/aws-hts/aws/");
source("scripts/R/lib/VCFs.R");
source("scripts/R/lib/Annotation.R");
source("scripts/R/lib/CNVs.R");

wdpath <- "/Users/seth/Work/nextGen_snp_analysis/deepseq_fdr_3_9/manu_peaks";
setwd(wdpath)
# rstart <- 33002251
# rend <- 33661647


rstart <- 33233467
rend <- 33365888

peakf <- read.table("FDs_3LpeakF_all.depth");
clips <- peakf[which(peakf$V2 == rstart):which(peakf$V2 == rend),]


depth <- peakf
window <- 10000
step <- 1000
rolldepth <-  rollapply(depth$V3,window,mean,by=step)
rollposn <-  rollapply(depth$V2,window,mean,by=step)
data<-as.data.frame(cbind(rollposn, rolldepth))
qplot(rollposn,rolldepth,data=data)
dev.print(png,file="3LPeakFlank_depthR_peakf.png",width=1200,height=800)



depth <- clips
window <- 10000
step <- 100
rolldepth <-  rollapply(depth$V3,window,mean,by=step)
rollposn <-  rollapply(depth$V2,window,mean,by=step)
data<-as.data.frame(cbind(rollposn, rolldepth))
qplot(rollposn,rolldepth,data=data)
dev.print(png,file="3LPeakFlank_depthR_clips.png",width=1200,height=800)


genes <- getGenesBM()
tmp <- genes[which(genes$chromosome_name == "3L"),]
tmp <- tmp[which(tmp$start_position >= rstart),]
tmp <- tmp[which(tmp$end_position <= rend),]
blockgenes <- cbind(tmp,clips=grepl("Clip",tmp$description))

vep <- getVep("./FDs_3LpeakF_all.vep","./FDs_3LpeakF_all.sev.vep")
tmp <- vep[which(vep$CHR == "3L"),]
tmp <- tmp[which(tmp$POS >= rstart),]
tmp <- tmp[which(tmp$POS <= rend),]
blockvep <- cbind(tmp,height=0)

blockvep <- blockvep[ which((!blockvep$Trans_consequence %in% c("UPSTREAM","DOWNSTREAM","INTRONIC","INTERGENIC","SYNONYMOUS_CODING"))),]
blockvep[ which((blockvep$Trans_consequence %in% c("NON_SYNONYMOUS_CODING","FRAMESHIFT_CODING","SPLICE_SITE","ESSENTIAL_SPLICE_SITE","STOP_GAINED","STOP_LOST"))),]$height=0
blockvep[ which((blockvep$Trans_consequence %in% c("FRAMESHIFT_CODING","SPLICE_SITE","ESSENTIAL_SPLICE_SITE","STOP_GAINED","STOP_LOST"))),]$height=10


geneLev <- 10
bheight <- 2


ggplot() +
  geom_line(data=data,aes(rollposn,rolldepth)) +
  geom_point(data=blockvep,aes(POS,geneLev*2+height,color=Max_consequence)) +
  geom_rect(data=blockgenes, aes(xmin=start_position, xmax=end_position,ymin=geneLev,ymax=(geneLev + (bheight * strand)),fill=clips)) +
  geom_text(data=blockgenes, aes(label=ensembl_gene_id, x = start_position+3000,y=(geneLev + ((bheight+2)*strand)),size=1,angle=(30*strand)))
dev.print(png,file="3LPeakFlank_depthR_clips_genes.png",width=1200,height=800)


candidates <- unique(blockvep[c("Uploaded_variation","POS","Gene","Max_consequence")])
candidates <- candidates[which(!candidates$Max_consequence %in% c("5PRIME_UTR","3PRIME_UTR","SPLICE_SITE","INTRONIC")),]
# candidates <- candidates[which(candidates$Gene %in% blockgenes[which(blockgenes$clips),"ensembl_gene_id"]),]
unique(merge(candidates,blockgenes[c("ensembl_gene_id","description")],by.x="Gene",by.y="ensembl_gene_id"))
candidates <- candidates[which(candidates$Max_consequence != "NON_SYNONYMOUS_CODING"),]
candidates <- unique(merge(candidates,blockgenes[c("ensembl_gene_id","description")],by.x="Gene",by.y="ensembl_gene_id"))
unique(merge(candidates,blockgenes[c("ensembl_gene_id","description")],by.x="Gene",by.y="ensembl_gene_id"))


expr_cand<- unique(blockvep[c("Uploaded_variation","POS","Gene","Max_consequence")])
expr_cand <- expr_cand[which(expr_cand$Gene %in% c("AGAP011789")),]
unique(merge(expr_cand,blockgenes[c("ensembl_gene_id","description")],by.x="Gene",by.y="ensembl_gene_id"))
