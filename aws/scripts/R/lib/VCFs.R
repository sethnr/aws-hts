## General variant-related annotation methods, parse VCFs, calc BAF, etc



## parse all VCFs in directory into tabix format
## vcfs <- getVcfs( "FD0.*.vcf.gz$")
getVcfs <- function(pattern) {
  vcfFiles <- list.files(path = ".", pattern = pattern)
  vcfTabs <- vector();
  vcfNames <- vector();
  vcfs <- vector();

  for (filename in vcfFiles) {
    write(filename, stderr());
    vcftab <- TabixFile(filename, index = paste(filename, "tbi", sep="."));
    vcfScan <- scanVcf(filename)
    
    write("unpacking",stderr())
    vcf <- unpackVcf(vcfScan, vcftab, info=TRUE, geno=FALSE)
    vcfs <- c(vcfs,vcf);
    write("vcfname",stderr())
    vcfNames <- c(vcfNames,gsub(".vcf.gz$","",filename));
    write("done",stderr());
  }
  rm(vcf, vcftab, vcfScan, vcfFiles, vcfTabs)
  
  names(vcfs) <- vcfNames
  return(vcfs)

}


#dbSNPalleles Assignment Chromosome Position Strand
# NOT YET WORKING!!! 
getInfo <- function(vcfs) {
  retDF <- data.frame()
  vcfNames <- names(vcfs)
  for(i in 1:length(vcfs)) {
    write(paste(i," ",vcfNames[[i]]),stderr())
    alleles <- paste(vcfs[[i]]$REF,vcfs[[i]]$ALT,sep="/");
    varName <- paste(vcfs[[i]]$CHROM,vcfs[[i]]$POS,alleles, sep = "_", collapse = NULL)
    retDF$NAME <- varName;
    retDF$dbSNPalleles <- alleles;
    retDF$Assignment <- alleles;
    retDF$Strand <- "+";
  }
  
  rm(varName, vcfNames, alleles);
  
  retDF[is.na(retDF)] <- 0
  colNames <- names(retDF)
  retDF<-cbind(t(as.data.frame(strsplit(as.character(retDF$Uploaded_variation),"[_]")))[,1:2])
  names(retDF) <- c("CHR","POS");
  retDF$POS <- as.numeric(as.character(retDF$POS))
  return(retDF)
}



## get B-allele frequencies from output of getVcf
## filter for minimum or maximum depth of coverage
## bafs <- getBafs(vcfs, minDP=3, maxDP=20)
getBafs <- function(vcfs, minDP=NULL, maxDP=NULL) {
  bafs <- data.frame()
  vcfNames <- names(vcfs)
  for(i in 1:length(vcfs)) {
    write(paste(i," ",names(vcfs)[[i]]),stderr())
    DP4l <- (length(vcfs[[i]]$INFO$DP4) / 4)
    dim(vcfs[[i]]$INFO$DP4) <- c(DP4l,4)
    
    vcfs[[i]]$INFO$BAF <- rowSums(vcfs[[i]]$INFO$DP4[,3:4],na.rm=FALSE)/vcfs[[i]]$INFO$DP
    varName <- paste(vcfs[[i]]$CHROM,vcfs[[i]]$POS,paste(vcfs[[i]]$REF,vcfs[[i]]$ALT,sep="/"), sep = "_", collapse = NULL)

    bafDF <- as.data.frame(list(varName, vcfs[[i]]$POS, vcfs[[i]]$INFO$BAF))
    names(bafDF) <- c("Uploaded_variation","POS",vcfNames[[i]]);
    #bafDF <- as.data.frame(list(varName, vcfs[[i]]$INFO$BAF))
    #names(bafDF) <- c("Uploaded_variation",vcfNames[[i]]);

    under <- data.frame(); over <- data.frame();
    if(exists("minDP") && (!is.null(minDP))) {
      pre <- length(bafDF[,1])
#      bafDF <- bafDF[which(vcfs[[i]]$INFO$DP >= minDP),];
      under <- bafDF[which(vcfs[[i]]$INFO$DP < minDP),"Uploaded_variation"]
#      write(paste("filtered: ",(length(bafDF[,1])/pre*100),"% >= ",minDP),stderr())      
      write(paste("filtered ",length(under)," < ",minDP),stderr())      
    }

    if(exists("maxDP") && (!is.null(maxDP))) {
      pre <- length(bafDF[,1])
#      bafDF <- bafDF[which(vcfs[[i]]$INFO$DP <= maxDP),];
      over <- bafDF[which(vcfs[[i]]$INFO$DP > maxDP),"Uploaded_variation"]
#      write(paste("filtered: ",(length(bafDF[,1])/pre*100),"% <= ",maxDP),stderr())      
      write(paste("filtered ",length(over)," > ",maxDP),stderr())      
    }

    pre <- length(bafDF[,1])      
    bafDF <- bafDF[(!bafDF$Uploaded_variation %in% under),];
    bafDF <- bafDF[(!bafDF$Uploaded_variation %in% over),];
    write(paste("filtered: ",pre," -> ",(length(bafDF[,1]))),stderr())      

#    if(length(bafs) > 0 ) { bafs <- merge(bafs,bafDF,all=TRUE,by="Uploaded_variation");}
    if(length(bafs) > 0 ) { bafs <- merge(bafs,bafDF,all=TRUE);}
    else { bafs <- bafDF}
  }
  
  rm(varName, DP4l, bafDF, vcfs);
  
  bafs[is.na(bafs)] <- 0
  return(bafs)
}


## parse all VCFs in directory into SNPmatrix
## vcfs <- getSM( "FD0.*.vcf.gz$")
getSMs <- function(pattern, grange) {
  vcfFiles <- list.files(path = ".", pattern = pattern)
  vcfTabs <- vector();
  vcfNames <- vector();
  SMs <- vector();

  for (filename in vcfFiles) {
    write(filename, stderr());
    vcftab <- TabixFile(filename, index = paste(filename, "tbi", sep="."));
    vcfSM <- vcf2sm(vcftab, grange, 7)
    
    SMs <- c(vcfs,vcf);
    write("vcfname",stderr())
    vcfNames <- c(vcfNames,gsub(".vcf.gz$","",filename));
    write("done",stderr());
  }
  rm(vcf, vcftab, vcfScan, vcfFiles, vcfTabs)
  
  names(SMs) <- vcfNames
  return(SMs)

}
