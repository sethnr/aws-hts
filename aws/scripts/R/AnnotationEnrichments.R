# go enrichment analysis, from biomart annotation sets from 'Annotation.R'

#genes2CC <- aggregate(x=genesGOCC$go_accession,by=list(genesGOCC$ensembl_transcript_id),FUN=function(x) {paste(x)})

genesGOCC <- genesGO[which(genesGO$namespace_1003 == "cellular_component"),]
genes2CC <- by(genesGOCC$go_accession,INDICES=list(genesGOCC$ensembl_transcript_id),FUN=function(x) {c(paste(x))})
# annFUN.gene2GO("CC",feasibleGenes = NULL,genes2CC)

genesGOBP <- genesGO[which(genesGO$namespace_1003 == "biological_process"),]
genes2BP <- by(genesGOBP$go_accession,INDICES=list(genesGOBP$ensembl_transcript_id),FUN=function(x) {c(paste(x))})

# annFUN.gene2GO("BP",feasibleGenes = NULL,genes2BP)

genesGOMF <- genesGO[which(genesGO$namespace_1003 == "molecular_function"),]
genes2MF <- by(genesGOMF$go_accession,INDICES=list(genesGOMF$ensembl_transcript_id),FUN=function(x) {c(paste(x))})

# annFUN.gene2GO("MF",feasibleGenes = NULL,genes2MF)

allGenes <- as.numeric(unique(genes$ensembl_transcript_id) %in% unique(vep$Feature))
names(allGenes) <- unique(genes$ensembl_transcript_id)


# GOFischer("CC",genes2CC)
GOFischer <- function(ontologyStr, genes2GO) {
  GOdata <- new("topGOdata", 
                  description = "Simple session", ontology = ontologyStr,
                  allGenes = allGenes, geneSel = function(x) {return (x==1)},
                  nodeSize = 10,
                  annot = annFUN.gene2GO, gene2GO=genes2GO)

  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  print(GenTable(GOdata, Fis = resultFisher, topNodes = 10))
}


GOFischer("CC",genes2CC)
GOFischer("MF",genes2MF)
GOFischer("BP",genes2BP)
