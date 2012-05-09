## General gene-related annotation methods, variant effect predictor, biomart, etc
library(biomaRt)


#############
#
# the following methods relate to biomart / biomaRt annotations
# get details of which mart using the following methods:
#    listMarts()
#    listDatasets(useMart("metazoa_mart_12"))
#
#############

## get data.frame of [attributes] for all genes (basis for following methods)
## getGenesBM(attributes) 
getGenesAtts <- function(attributes) {
  ensembl <- useMart("metazoa_mart_13", "agambiae_eg_gene")
  genes <- getBM(attributes=attributes, mart = ensembl)
  return(genes)
}


## get data.frame of gene / desc / location info for all genes
## getGenesBM()
getGenesBM <- function() {
  attributes  <- c("ensembl_gene_id",
                  "ensembl_transcript_id",
                  "description",
                  "chromosome_name",
                  "start_position",
                  "end_position",
                  "strand")

  return(getGenesAtts(attributes))
}

## get data.frame of gene / go_term / desc info for all genes
## getGenesBM()
getGenesGO <- function() {
  attributes <- c("ensembl_gene_id",
                  "ensembl_transcript_id",
                  "go_accession",
                  "namespace_1003",
                  "name_1006")

  return(getGenesAtts(attributes))
}

## get data.frame of gene / interpro_id / desc info for all genes
## getGenesIpr()
getGenesIpr <- function() {
  attributes <- c("ensembl_gene_id",
                  "ensembl_transcript_id",
                  "interpro",
                  "interpro_description"
                  )
  return(getGenesAtts(attributes))
}

#############
#
#    the following methods relate to Variant Effect Predictor annotations
#    see website for further details:
#       http://www.ensembl.org/info/docs/variation/vep/index.html
#
#############

## parse output of VEP
## acts on up to two files, first called with normal sub, second with --severity to get max_consequence.
## vep <- getVep("./split/all.vep.txt",maxSevFile="./split/all.vep.sev.txt")
getVep <- function(consFile, maxSevFile=NULL) {
  write(paste("getVep(",consFile," ",maxSevFile,")"),stderr())

  write(paste("reading: ",consFile),stderr())
  vep <- read.table(consFile,sep="\t")
  vep<-cbind(t(as.data.frame(strsplit(as.character(vep$V1),"[_]")))[,1:2],vep)
  names(vep) <- c("CHR","POS","Uploaded_variation", "Location","Allele",  "Gene", "Feature", "Feature_type",    "Trans_consequence",     "cDNA_position",   "CDS_position",    "Protein_position",        "Amino_acids",     "Codons",  "Existing_variation",      "Extra")
  vep$POS <- as.factor(vep$POS)

  

  if(exists("maxSevFile") && (!is.null(maxSevFile))) {
    write(paste("reading: ",maxSevFile),stderr())
    max <- read.table(maxSevFile,sep="\t")
    max <- data.frame(max[1],max[7])
    names(max) <- c("Uploaded_variation","Max_consequence")
    write(paste("merging: ",consFile," + ",maxSevFile),stderr())
    vep <- merge(vep,max,by="Uploaded_variation",all=TRUE)
    rm(max)
  }
  vep$POS <- as.numeric(as.character(vep$POS))
  return(vep)
} 


# calculate DN:DS (non synonymous:synonymous ratios across a VEP file
# Nb, acts on output of getVcf with maxSevFile (see getVep(consFile, maxSevFile))
#
DNDS <- function(data) 
{
  nsCons <- c("NON_SYNONYMOUS_CODING","FRAMESHIFT_CODING","SPLICE_SITE","ESSENTIAL_SPLICE_SITE","STOP_GAINED","STOP_LOST")
  sCons <- c("SYNONYMOUS_CODING","3PRIME_UTR","5PRIME_UTR","WITHIN_NON_CODING_GENE","CODING_UNKNOWN")
  ncCons <- c("INTRONIC","UPSTREAM","DOWNSTREAM","INTERGENIC")
  allCons <- names(summary(data$Max_consequence))
  error <- allCons[!(allCons %in% c(nsCons,sCons,ncCons))]
  write(paste("warning:", error, "not tested"),stderr())
  
  cols <- c("Uploaded_variation","Feature","Max_consequence")
  data2 <- data[ which(!(data$Max_consequence %in% ncCons)),]
  
  nonsyn <- unique(data2[ which((data2$Max_consequence %in% nsCons)),cols])
  syn <- unique(data2[ which((data2$Max_consequence %in% sCons)),cols])
  
  summary(data2); summary(syn); summary(nonsyn)

  write("aggregating across genes",stderr())
  
  synL <- aggregate(syn$Uploaded_variation, by=list(syn$Feature), FUN=function(x) {length(unique(x))})
  names(synL) <- c("Feature","n.S");
  nonSynL <- aggregate(nonsyn$Uploaded_variation, by=list(nonsyn$Feature), FUN=function(x) {length(unique(x))})
  names(nonSynL) <- c("Feature","n.NS");

  write(paste("synL",dim(synL),"\n","nonSynL",dim(nonSynL)),stderr())

  dnds <- merge(nonSynL,synL,all=TRUE)
  dnds <- cbind(dnds,dnds = (dnds$n.NS/dnds$n.S))
  return(dnds)
}

