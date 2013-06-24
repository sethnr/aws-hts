library(topGO)
library(reshape2)
library(biomaRt)


options(width = 120)

splitGOs <- function() {
  agamGO <- read.delim("all_dmel_GO.txt", sep="\t",header=F)
  agamGO[[2]] <- as.character(agamGO[[2]])
  head(agamGO)
  GOno <- agamGO[which(agamGO[[3]]==""),1:2]
  GOcc <- agamGO[which(agamGO[[3]] == "cellular_component"),1:2]
  GOmf <- agamGO[which(agamGO[[3]] == "molecular_function"),1:2]
  GObp <- agamGO[which(agamGO[[3]] == "biological_process"),1:2]
  GOcc <- rbind(GOno,GOcc)
  GOmf <- rbind(GOno,GOmf)
  GObp <- rbind(GOno,GObp)
  GOcc <- as.data.frame(aggregate(GOcc[[2]], by=list(GOcc[[1]]),collapse=';',FUN=paste))
  GOmf <- as.data.frame(aggregate(GOmf[[2]], by=list(GOmf[[1]]),collapse=';',FUN=paste))
  GObp <- as.data.frame(aggregate(GObp[[2]], by=list(GObp[[1]]),collapse=';',FUN=paste))
  write.table(GOcc, "all_GO.cc.txt", row.names=F, col.names=F,sep="\t",quote=F)
  write.table(GOmf, "all_GO.mf.txt", row.names=F, col.names=F,sep="\t",quote=F)
  write.table(GObp, "all_GO.bp.txt", row.names=F, col.names=F,sep="\t",quote=F)
}

getBiomarts <- function(mart_name) {
 listMarts()
 ensembl=useMart("ensembl")
 #... finish writing this shit when there's a working connection at pasteur. 
 #Proxy bollocks again...
  
}



runTopGO <- function(ontology, gene2GO, allGenes, cluster) {
  inCluster <- as.numeric(allGenes[[1]] %in% cluster[[1]])
  names(inCluster) <- allGenes[[1]]
  f1 <- function(x) {return(x==1)}
  GOdataCC1 <- new("topGOdata", 
                     description = "dmel topGO analysis", ontology = ontology,
                     allGenes = inCluster, geneSel = f1,
                     nodeSize = 10,
                     annot = annFUN.gene2GO, gene2GO=gene2GO)


  resultFisherCC1 <- runTest(GOdataCC1, algorithm = "classic", statistic = "fisher")

  print(paste("GO results cluster ",ontology),outfile)
  table <- GenTable(GOdataCC1, Fis = resultFisherCC1, topNodes = 10)
  write.table(as.data.frame(table),sep="\t", quote=F, row.names=F, outfile)
  
  sigGOs <- table[which(table$Fis <= 0.05),1]
  write(paste("sig GO terms: ",sigGOs),stderr())
  GOM <- melt(gene2GO)
  names(GOM) <- c("GO","GENE")
  GOgenes <- GOM[which(GOM$GO %in% sigGOs),2]
  print(paste("genes with sig. GO terms ",length(GOgenes)),outfile)
#  GOgenes <- GOM[which(GOM$GENE %in% inCluster),2]
#  print(paste("in cluster ",dim(GOgenes)),outfile)
  
  print(paste("genes with sig. GO terms ",ontology),outfile)
  write.table(cluster[which(cluster[[1]] %in% GOgenes),],sep="\t", quote=F, row.names=F, outfile)
  
}


splitGOs()

GOCC <- readMappings("all_GO.cc.txt",IDsep=";")
GOBP <- readMappings("all_GO.bp.txt",IDsep=";")
GOMF <- readMappings("all_GO.mf.txt",IDsep=";")


cmd_args = commandArgs(trailingOnly = T);
allGenesFile <- "./allGenes.txt"
outfile <- stdout()
for (arg in cmd_args) {
  write(arg,stderr());
  arg <- strsplit(arg,"=")[[1]];
  if(arg[1] == "allGenes") allGenesFile <- as.character(arg[2]);
  if(arg[1] == "genes") allGenesFile <- as.character(arg[2]);
  if(arg[1] == "out") outfile <- file(as.character(arg[2]), "w");
  if(arg[1] == "out") open(outfile);

}

write(paste("reading allGenes: ",allGenesFile),stderr())
allGenes <- as.data.frame(read.table(allGenesFile,header=T,sep="\t"))

tmpFileName <- "/tmp/runTopGO.TMP.txt"

write("reading cluster",stderr())
tmpFile <- file(tmpFileName, "w")
fin <- file("stdin")
open(fin)
while(length(line <- readLines(fin,n=1)) > 0) {
  #nb, scanVcf will flip if given trailing tabs, chomp everything. 
  write(line, tmpFile)
}           
close(fin)
close(tmpFile)
write("finished writing tmp cluster",stderr())

clusterTab <- read.table(tmpFileName, header=T,sep="\t")
#write(head(clusterTab),stderr())
write("reading cluster",stderr())
cluster <- as.data.frame(clusterTab)
#write(head(cluster),stderr())
write("finished reading cluster",stderr())

runTopGO("CC",GOCC,allGenes,cluster) 
runTopGO("MF",GOMF,allGenes,cluster)
runTopGO("BP",GOBP,allGenes,cluster)
if (outfile != stdout()) close(outfile)
