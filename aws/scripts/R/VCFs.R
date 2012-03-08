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
    bafDF <- as.data.frame(list(varName, vcfs[[i]]$INFO$BAF))
    names(bafDF) <- c("Uploaded_variation",vcfNames[[i]]);

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

    if(length(bafs) > 0 ) { bafs <- merge(bafs,bafDF,all=TRUE,by="Uploaded_variation");}
    else { bafs <- bafDF}
  }
  
  rm(varName, DP4l, bafDF, vcfs);
  
  bafs[is.na(bafs)] <- 0
  bafs <- cbind(bafs,Contrast=abs(bafs$FD03_high.peaks - bafs$FD03_zero.peaks))
  return(bafs)
}

