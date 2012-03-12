## Annotation methods relating to CNV, raw depth, etc
# to add:
#  normalisation of raw results
#  segmentation analysis via ultrasome

getRawDP <- function(vcfs) {
  dps <- data.frame()
  vcfNames <- names(vcfs)
  for(i in 1:length(vcfs)) {
    write(paste(i," ",vcfNames[[i]]),stderr())
    varName <- paste(vcfs[[i]]$CHROM,vcfs[[i]]$POS,paste(vcfs[[i]]$REF,vcfs[[i]]$ALT,sep="/"), sep = "_", collapse = NULL)
    dpDF <- as.data.frame(list(varName, vcfs[[i]]$INFO$DP))
    names(dpDF) <- c("Uploaded_variation",vcfNames[[i]]);
    if(length(dps) > 0 ) { dps <- merge(dps,dpDF,all=TRUE,by="Uploaded_variation");}
    else { dps <- dpDF}
  }
  
  rm(dpDF, varName, vcfNames);
  
  dps[is.na(dps)] <- 0
  colNames <- names(dps)
  dps<-cbind(t(as.data.frame(strsplit(as.character(dps$Uploaded_variation),"[_]")))[,1:2],dps)
  names(dps) <- c("CHR","POS",colNames);
  dps$POS <- as.numeric(as.character(dps$POS))
  return(dps)
}


# parse out depths from VCF file, run summary(depth) and return result
# 
getDPStats <- function(vcfs) {
  dp <- vector();
  for(i in 1:length(vcfs)) {
    dp <- c(dp,vcfs[[i]]$INFO$DP);
  }
  summary(dp)
}

