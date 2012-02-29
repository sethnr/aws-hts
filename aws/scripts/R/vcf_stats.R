# use ggtools to read in VCF, calc fst across csome.
# having previously run RSetUp.R

library(GGtools)
library(Rsamtools)

wdpath <- "/Users/seth/Work/nextGen_snp_analysis/deepseq_fdr_3_9/";
setwd(wdpath)
vcfFiles <- list.files(path = ".", pattern = "FD0*.vcf.gz$")
vcfTabs <- vector();
vcfNames <- vector();
vcfs <- vector();

myGR <- gr2L;

for (filename in vcfFiles) {
#  vcftab <- TabixFile(filename, index = paste(filename, "tbi", sep="."));
#  vcf <- scanVcf(vcfTabs, param=myGR)
str(filename);
  vcftab <- TabixFile(filename, index = paste(filename, "tbi", sep="."));
  vcfScan <- scanVcf(filename)

str("unpacking")
  vcf <- unpackVcf(vcfScan, vcftab, info=TRUE, geno=FALSE)
  vcfs <- c(vcfs,vcf);
str("vcfname")
  vcfNames <- c(vcfNames,gsub(".vcf.gz$","",filename));
str("done");
}

names(vcfs) <- vcfNames

bafs <- data.frame()
for(i in 1:length(vcfs)) {
  DP4l <- (length(vcfs[[i]]$INFO$DP4) / 4)
  dim(vcfs[[i]]$INFO$DP4) <- c(DP4l,4)
  vcfs[[i]]$INFO$BAF <- rowSums(vcfs[[i]]$INFO$DP4[,3:4],na.rm=FALSE)/vcfs[[i]]$INFO$DP
  varName <- paste(vcfs[[i]]$CHROM,vcfs[[i]]$POS,paste(vcfs[[i]]$REF,vcfs[[i]]$ALT,sep="/"), sep = "_", collapse = NULL)
  bafDF <- as.data.frame(list(varName, vcfs[[i]]$INFO$BAF))
  names(bafDF) <- c("Uploaded_variation",vcfNames[[i]]);
  if(length(bafs) > 0 ) { bafs <- merge(bafs,bafDF,all=TRUE,by="Uploaded_variation");}
  else { bafs <- bafDF}
}
bafs[is.na(bafs)] <- 0
bafs <- cbind(bafs,Contrast=abs(bafs$FD03_high.peaks - bafs$FD03_zero.peaks))



#plot rolling windows of 100snps, 5 snps step size
library(zoo)
#plot(x=rollapply(bafs$POS,100,max, by=5), y=rollapply(abs(bafs$FD03_high.peaks - bafs$FD03_zero.peaks),100,mean, by=5), xlim=c(33002251,33661647))


vep <- read.table("./split/all.vep.sev.txt",sep="\t")
# head(t(strsplit(as.character(vep$V2),":",fixed-TRUE)))
# cons<-cbind(t(as.data.frame(strsplit(as.character(vep$V1),"[_]")))[,1:2],as.data.frame(vep$V7))
cons <- data.frame(vep[1],vep[7])
names(cons) <- c("Uploaded_variation","Max_consequence")
# names(cons) <- c("CHR","POS","Max_consequence")

vep <- read.table("./split/all.vep.txt",sep="\t")
vep<-cbind(t(as.data.frame(strsplit(as.character(vep$V1),"[_]")))[,1:2],vep)
names(vep) <- c("CHR","POS","Uploaded_variation", "Location","Allele",  "Gene", "Feature", "Feature_type",    "Trans_consequence",     "cDNA_position",   "CDS_position",    "Protein_position",        "Amino_acids",     "Codons",  "Existing_variation",      "Extra")
vep$POS <- as.factor(vep$POS)
vep <- merge(vep,cons,by="Uploaded_variation",all=TRUE)


data <- merge(vep, bafs, all.y=TRUE ,by="Uploaded_variation")
#change factor to numeric values to numeric values
data$POS <- as.numeric(as.character(data$POS))
rm(vep, cons)



coding_data <- data[ which(!(data$Max_consequence %in% c("UPSTREAM","DOWNSTREAM","INTRONIC"))),]
#remove any silent mutations:
nonsyn_data <- data[ which((data$Max_consequence %in% c("NON_SYNONYMOUS_CODING","FRAMESHIFT_CODING","SPLICE_SITE","ESSENTIAL_SPLICE_SITE","STOP_GAINED","STOP_LOST"))),]

# bind contrast to DF
#coding_data <- cbind(coding_data,Contrast=abs(coding_data$FD03_high.peaks - coding_data$FD03_zero.peaks))
qplot(POS,y=Contrast,data=coding_data,xlim=c(33002251,33661647),color=Max_consequence)



# listMarts()
# listDatasets(useMart("metazoa_mart_12"))
getGenesBM <- function() {
  ensembl <- useMart("metazoa_mart_12", "agambiae_eg_gene")

  attributes <- c("ensembl_gene_id",
                  "ensembl_transcript_id",
                  "description",
                  "chromosome_name",
                  "start_position",
                  "end_position",
                  "strand",
                  "interpro",
                  "interpro_description",
                  "go_accession",
                  "go_domain",
                  "name_1006")
  
  genes <- getBM(attributes=attributes, mart = ensembl)
  return(genes)
}
  go <- unique(as.data.frame(cbind(genes$ensembl_transcript_id,genes$go_accession,genes$name_1006)))
  ipr <- unique(as.data.frame(cbind(genes$ensembl_transcript_id,genes$interpro,genes$interpro_description)))
  

DNDS <- function (data) 
{
  nsCons <- c("NON_SYNONYMOUS_CODING","FRAMESHIFT_CODING","SPLICE_SITE","ESSENTIAL_SPLICE_SITE","STOP_GAINED","STOP_LOST")
  sCons <- c("SYNONYMOUS_CODING","3PRIME_UTR","5PRIME_UTR","WITHIN_NON_CODING_GENE","CODING_UNKNOWN")
  ncCons <- c("INTRONIC","UPSTREAM","DOWNSTREAM","INTERGENIC")
  allCons <- names(summary(data$Max_consequence))
  error <- allCons[!(allCons %in% c(nsCons,sCons,ncCons))]
  write(paste("warning:", error, "not tested"),stderr())
  
  
                                        # cols <- c("Uploaded_variation","Feature","Trans_consequence","Max_consequence")
  cols <- c("Uploaded_variation","Feature","Max_consequence")
  data2 <- data[ which(!(data$Max_consequence %in% ncCons)),]
  
  nonsyn <- unique(data2[ which((data2$Max_consequence %in% nsCons)),cols])
  syn <- unique(data2[ which((data2$Max_consequence %in% sCons)),cols])
  
  summary(data2); summary(syn); summary(nonsyn)

  write("aggregating across genes",stderr())
 
  #aggregate(syn$Uploaded_variation, by=list(syn$Feature), FUN=function(x) {length(unique(x))})

#  syns <- aggregate(syn$Uploaded_variation, by=list(syn$Feature), FUN=function(x) {paste(unique(x))})
#  names(syns) <- c("Feature","S");
#  nonSyns <- aggregate(nonsyn$Uploaded_variation, by=list(nonsyn$Feature), FUN=function(x) {paste(unique(x))})
#  names(nonSyns) <- c("Feature","NS");
  
  synL <- aggregate(syn$Uploaded_variation, by=list(syn$Feature), FUN=function(x) {length(unique(x))})
  names(synL) <- c("Feature","n.S");
  nonSynL <- aggregate(nonsyn$Uploaded_variation, by=list(nonsyn$Feature), FUN=function(x) {length(unique(x))})
  names(nonSynL) <- c("Feature","n.NS");

  write(paste("synL",dim(synL),"\n","nonSynL",dim(nonSynL)),stderr())

  
#  dnds <- merge(merge(nonSynL,nonSyns),merge(synL,syns))
  dnds <- merge(nonSynL,synL,all=TRUE)
  dnds <- cbind(dnds,dnds = (dnds$n.NS/dnds$n.S))
  return(dnds)
}
