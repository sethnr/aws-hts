



genes <- getGenesBM()
genesGO <- getGenesGO()
genesIpr <- getGenesIpr()

# get non-synonymous vars & merge with gene GO terms
diffNonsyn <- merge(nonsyn_data[which(nonsyn_data$Contrast >= 0.6),c("Uploaded_variation","CHR","POS","Gene","Feature","Trans_consequence","Max_consequence","Contrast")], genesGO, all.x=TRUE, by.x="Feature", by.y="ensembl_transcript_id")





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


