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

#getVCF <- function(filename, params) {
getVCF <- function(filename) {
  write(filename, stderr());
#  write(params, stderr());
  vcftab <- TabixFile(filename, index = paste(filename, "tbi", sep="."));
#  vcfScan <- scanVcf(filename, params=params)
  vcfScan <- scanVcf(filename)
  
  write("unpacking",stderr())
  vcf <- unpackVcf(vcfScan, vcftab, info=TRUE, geno=TRUE)
  write("vcfname",stderr())
  write("done",stderr());  
  vcf
}

getBafsMultiFile <- function(vcf) {
  noSamples <- dim(vcf[[1]]$GENO$DP4)[2];
  noVars <- dim(vcf[[1]]$GENO$DP4)[1];
  bafs <- character()
  for(s in 1:noSamples) {
    dp4 <- vcf[[1]]$GENO$DP4[,s,];
    baf <- apply(dp4,1,function(x) {sum(x[3:4]) / sum(x);});
    if(length(bafs) > 0) {write("binding",stderr())
                          bafs<-cbind(bafs,baf);}
    else {bafs <- as.matrix(baf);}

    write(dim(bafs),stderr())
    write(bafs[1:3,],stderr())
    
  }
  vcf[[1]]$GENO$BAF <- bafs;
  vcf
}


calcThetaVcf <- function(vcf, start=NULL, end=NULL, chr=NULL) {
  noSamples <- dim(vcf[[1]]$GENO$DP4)[2];
  thetas <- vector()
  bafs <- vcf[[1]]$GENO$BAF;
  bafs <- data.frame(vcf[[1]]$CHR, vcf[[1]]$POS, vcf[[1]]$GENO$BAF)
  dps <- vcf[[1]]$GENO$DP;
  dps <- data.frame(vcf[[1]]$CHR, vcf[[1]]$POS, vcf[[1]]$GENO$DP)
  names(bafs)[1:2] <- c("CHR","POS");
  names(dps)[1:2] <- c("CHR","POS");
  
  if(!is.null(start) || !is.null(end)) {
    if (exists("chr")) { 
        bafs <- bafs[which(bafs$CHR == chr),];
        dps <- dps[which(dps$CHR == chr),];
    }  
    bafs <- bafs[which(bafs$POS <= end),];
    bafs <- bafs[which(bafs$POS >= start),];
    dps <- dps[which(dps$POS <= end),];
    dps <- dps[which(dps$POS >= start),];
  }
  sites <- max(bafs$POS) - min(bafs$POS);

  #values for total theta: 
  sumMeanDepth = 0; sumSeg = 0;

  for(s in 3:(noSamples+2)) {
    dp <- dps[,s];
    dp <- dp[which(!is.na(dp))]
    meanDp <- mean(dp);
    baf <- bafs[,s];
    segSites <- length(baf[which(baf > 0)])
  
    sumMeanDepth <- sumMeanDepth + meanDp
    
    theta <- segSites  / nthHarmonic((meanDp -1))
    thetas <- c(thetas, theta)
    }  

  theta <- dim(bafs)[[1]]  / nthHarmonic((sumMeanDepth -1))
  thetas <- c(thetas, theta)

  thetas
}
nthHarmonic <- function(n) { h<-0; for (i in 1:n) { h<-(h+(1/i));}; h; }




## get B-allele frequencies from output of getVcf
## filter for minimum or maximum depth of coverage
## bafs <- getBafs(vcfs, minDP=3, maxDP=20)
getBafs <- function(vcf, minDP=NULL, maxDP=NULL) {
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
getSM <- function(filename, grange) {
  write(filename, stderr());
  #  write(params, stderr());
  vcftab <- TabixFile(filename, index = paste(filename, "tbi", sep="."));
  vcfSM <- vcf2sm(vcftab, gr=grange, nmetacol=9L)
      
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
