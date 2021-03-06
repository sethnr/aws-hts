# load required libraries
#library("zoo")
library("ShortRead")
library("snpStats");
library("GGtools");
library("VariantAnnotation")
library("GenomicRanges")

########
# standard modules (shift to dedicated library)
########


getVCF <- function(file_name) {
  write(file_name, stderr());
#  write(params, stderr());
  
  params <- ScanVcfParam(trimEmpty=T)
  vcftab <- TabixFile(file_name, index = paste(file_name, "tbi", sep="."));
  vcftab
  # vcfScan <- scanVcf(filename, params=params)
  write("scanning",stderr())

  vcfScan <- scanVcf(file_name)
  
  write("unpacking",stderr())
  vcf <- unpackVcf(vcfScan, vcftab, info=TRUE, geno=TRUE)
  write("vcfname",stderr())
  write("done",stderr());  
  vcf
}

getBafsMultiFile <- function(vcf) {
  
  for (i in 1:length(vcf)) {
#    write(paste("getting BAFs ",i),stderr())
    noSamples <- dim(vcf[[i]]$GENO$DP4)[2];
    noVars <- dim(vcf[[i]]$GENO$DP4)[1];
    bafs <- character()
    for(s in 1:noSamples) {
      dp4 <- vcf[[i]]$GENO$DP4[,s,];
      baf <- apply(dp4,1,function(x) {sum(x[3:4]) / sum(x);});
      if(length(bafs) > 0) { #write("binding",stderr())
                            bafs<-cbind(bafs,baf);}
      else {bafs <- as.matrix(baf);}
      
#      write(dim(bafs),stderr())
#      write(bafs[1:3,],stderr())
      
    }
    colnames(bafs) <- colnames(vcf[[i]]$GENO$DP4)
    vcf[[i]]$GENO$BAF <- bafs;
  }
 vcf
}

getUgroups <- function(colnames) {
  groups <- colnames
  groups <- matrix(unlist(strsplit(as.character(groups),split="_")),nrow=2)[1,]
  
  #CODE TO HANDLE MULTIPLE GROUPS (i.e. NAs)
  ugroups <- list()
  ugroups$vals <- list()
  ugroups$groups <- unique(groups)
  for (g in 1:length(ugroups)) {
    ugroups$vals[g] <- list(which(groups == ugroups$groups[g]))
  }
  ugroups
}


getMaxReadShift <- function(vcf) {  
  for (i in 1:length(vcf)) {
#    write(paste("getting max read shift ",i),stderr())
    
    noSamples <- dim(vcf[[i]]$GENO$DP4)[2];
    noVars <- dim(vcf[[i]]$GENO$DP4)[1];
    #  MBCs <- vector()
#    groups <- colnames(vcf[[i]]$GENO$BAF)
#    groups <- matrix(unlist(strsplit(as.character(groups),split="_")),nrow=2)[1,]
    
    #CODE TO HANDLE MULTIPLE GROUPS (i.e. NAs)
#    ugroups <- list()
#    ugroups$vals <- list()
#    ugroups$groups <- unique(groups)
#    for (g in 1:length(ugroups)) {
#      ugroups$vals[g] <- list(which(groups == ugroups$groups[g]))
#    }
    ugroups <- getUgroups(colnames(vcf[[i]]$GENO$BAF));

    #  MBCs <- vector()
    for (g in 1:length(ugroups$groups)) {  
      #    maxBafContrast <- apply(bafs,1,function(x) {max(x[which(!is.na(x))]) - min(x[which(!is.na(x))])});
      groupCols <- c(ugroups$vals[[g]])
      
      groupTotalBAF <- apply(vcf[[i]]$GENO$BAF[,groupCols],1,sum)/3
#      write(length(groupTotalBAF),stderr())
      inv <- which(groupTotalBAF > 0.5)
      groupMafs <- vcf[[i]]$GENO$BAF[,groupCols];
      groupMafs[inv,] <- 1-groupMafs[inv,]
      
      groupTotalDP <- apply(vcf[[i]]$GENO$DP[,groupCols],1,sum)
      groupDPs <- vcf[[i]]$GENO$DP[,groupCols]

      #minor allele depth
      groupMaDP <- groupMafs * groupDPs
#      write(dim(groupMaDP),stderr())
      maxReadShift <- apply(groupMaDP,1,function(x) {max(x) - min(x)});
      
      
      if(exists("MRSs")) {
        MRSs <- cbind(MRSs,maxReadShift)
      }
      else {
        MRSs <- maxReadShift
      }
    }
    #colnames(MRSs) <- ugroups$groups
    vcf[[i]]$INFO$MRS<- MRSs;
    rm(MRSs)
  }
  vcf
}



getMaxBafContrast <- function(vcf) {
  
  for (i in 1:length(vcf)) {
#    write(paste("getting max contrasts ",i),stderr())
    
    noSamples <- dim(vcf[[i]]$GENO$DP4)[2];
    noVars <- dim(vcf[[i]]$GENO$DP4)[1];
    #  MBCs <- vector()
    groups <- colnames(vcf[[i]]$GENO$BAF)
    groups <- matrix(unlist(strsplit(as.character(groups),split="_")),nrow=2)[1,]
    
    #CODE TO HANDLE MULTIPLE GROUPS (i.e. NAs)
    ugroups <- list()
    ugroups$vals <- list()
    ugroups$groups <- unique(groups)
    for (g in 1:length(ugroups)) {
      ugroups$vals[g] <- list(which(groups == ugroups$groups[g]))
    }
    
    #  MBCs <- vector()
    for (g in 1:length(ugroups$groups)) {  
      bafs <- vcf[[i]]$GENO$BAF;
      #    maxBafContrast <- apply(bafs,1,function(x) {max(x[which(!is.na(x))]) - min(x[which(!is.na(x))])});
      maxBafContrast <- apply(bafs,1,function(x) {max(x[c(ugroups$vals[[g]])]) - min(x[c(ugroups$vals[[g]])])});
      
      if(exists("MBCs")) {
        MBCs <- cbind(MBCs,maxBafContrast)
      }
      else {
        MBCs <- maxBafContrast
      }
    }
    #colnames(MBCs) <- ugroups$groups
    vcf[[i]]$INFO$MC <- MBCs;
    rm(MBCs)
  }
  vcf
}

getNonsenseRows <- function(x, nonsense_cons) {
  nonsense <- FALSE
  for (n in 1:length(nonsense_cons)) {
    for (c in 1:length(x[[1]])) {
      if (x[[1]][c] == nonsense_cons[n]) {nonsense <- TRUE}
    }
  }
  nonsense
}


printValidVars <- function(vcf, minMC, minRS, CSQs, genes) {
  write(paste(CSQs),stderr())
  for (i in 1:length(vcf)) {
    # write(paste("printing valids ",i),stderr())
    csqs <- vcf[[i]]$INFO$CSQ
    good_csq <- sapply(csqs,FUN = getNonsenseRows, nonsense_cons=CSQs)
    csqs <- sapply(csqs,FUN=function(x) {paste(x,sep=":",collapse=":")})

#    if(exists("genes")) {
#      vars <- cbind(vcf[[i]]$CHROM ,
#                    vcf[[i]]$POS ,
#                    vcf[[i]]$ALT)
#      colnames(vars) <- c("CHR","POS","ALLELE")
#      vars <- merge(vars,genes)      
#      }
    
    mcs <- vcf[[i]]$INFO$MC;
    mcs[which(is.na(mcs))] <- 0
    if (is.matrix(mcs)) {mc <- apply(mcs,1,FUN=max)}
    else {mc <- mcs}
    
    good_mc <- mc >= minMC
    
    rss <- vcf[[i]]$INFO$MRS;
    rss[which(is.na(rss))] <- 0
    if (is.matrix(rss)) {rs <- apply(rss,1,FUN=max)}
    else {rs <- rss}
#    rs <- apply(rss,1,FUN=max)
    good_mrs <- rs >= minRS
    
    good <- apply(cbind(good_mc, good_mrs, good_csq), 1, FUN=function(x) {x[1] && x[2] && x[3]})
    good <- which(good)
    
    out <- cbind(
      vcf[[i]]$CHROM ,
      vcf[[i]]$POS ,
      vcf[[i]]$REF ,
      vcf[[i]]$ALT ,
      vcf[[i]]$INFO$MRS ,
      vcf[[i]]$INFO$MC ,
      vcf[[i]]$INFO$FS ,
      csqs
    )
#    names <- c("CHR","POS","REF","ALT",colnames(vcf[[i]]$INFO$MRS),colnames(vcf[[i]]$INFO$MC),"CSQ")
    out <- matrix(out[good,],ncol=dim(out)[[2]])
#    colnames(out) <- names
#    if (i==1) {
#      write.table(out,col.names=T, row.names=F, quote=F, sep="\t",stdout())
#    }
#    else {
      write.table(out,col.names=F, row.names=F, quote=F, sep="\t",stdout())
#      }
     }  
}




calcThetaBAFs <- function(x, w_start, w_end) {
  # x = BAF/DP matrix
  if(!is.matrix(x)) {write("I was expecting a matrix", stderr());
                    }

  noSamples <- (dim(x)[2]-2) /2;
  noSamples <- (dim(x)[2]-2) /2;
  thetas <- vector()
#  thetas <- c(thetas, )
  posns <- as.numeric(x[,2])
#  thetas <- c(thetas,x[1,1],min(m[,2]), max(m[,2]))  
  if(exists("w_start") && exists("w_end")) {
    thetas <- c(thetas,x[1,1],w_start,w_end)  
  }
  else {
    thetas <- c(thetas,x[1,1],min(posns), max(posns))  
  }
  #values for total theta: 
  sumMeanDepth = 0; sumSeg = 0;
  
  for(s in 1:noSamples) {
    dp <- x[,s+2+noSamples];
    meanDp <- mean(as.numeric(dp),na.rm=T);
      # limit  depth to max pool size
    sumMeanDepth <- sumMeanDepth + meanDp

    baf <- x[,s+2];
    baf[which(is.na(baf))] <- 0 
    segSites <- length(baf[which(baf > 0)])
    
    theta <- segSites  / nthHarmonic((meanDp -1))
#    print(paste(segSites, "/ nthHarmonic(", meanDp ,") -> ",theta),stderr());
    thetas <- c(thetas, theta)
    rm(dp,meanDp,baf,segSites,theta)
    gc()
    }  
#  print(x[which(is.na(x))],stderr())
#  print(is.matrix(x),stderr())
#  print(is.data.frame(x),stderr())
  x[which(is.na(x))] <- 0 
#  print(x[which(is.na(x))],stderr())

  bafTotal <- apply(x[,3:(noSamples+2)],1,function(x) {sum(as.numeric(x))})
#  print(head(bafTotal),stderr())
  segTotal <- length(bafTotal[which(bafTotal > 0)])
#  print(head(segTotal),stderr())
  
  theta <- segTotal  / nthHarmonic((sumMeanDepth -1))
  thetas <- c(thetas, theta)
  rm(bafTotal, segTotal, theta)  
  write(cat(thetas,sep="\t"),stdout())
  rm(thetas)
}



chomp <- function (x) sub("\\s+$", "", x)







coding_e <- c("SYNONYMOUS_CODING","SPLICE_SITE","5_PRIME_UTR","3_PRIME_UTR",
                "TRANSCRIPT_ABLATION","ESSENTIAL_SPLICE_SITE","STOP_GAINED",
                "FRAMESHIFT_CODING","STOP_LOST","NON_SYNONYMOUS_CODING",
                "REGULATORY_REGION","REGULATORY_REGION_ABLATION","REGULATORY_REGION_AMPLIFICATION")
missense_e <- c("TRANSCRIPT_ABLATION","ESSENTIAL_SPLICE_SITE","STOP_GAINED",
                "FRAMESHIFT_CODING","STOP_LOST","NON_SYNONYMOUS_CODING",
                "REGULATORY_REGION","REGULATORY_REGION_ABLATION","REGULATORY_REGION_AMPLIFICATION")
nonsense_e <- c("TRANSCRIPT_ABLATION","ESSENTIAL_SPLICE_SITE","STOP_GAINED",
                "FRAMESHIFT_CODING","STOP_LOST")



# parse arguments from command line
cmd_args = commandArgs(trailingOnly = T);

#windowsFile <- "Fd03_a.theta.sd.Q3_3"
#vcfFileName <- "Fd.combined.vep.vcf"
#setwd("~/Scratch/Fd_VEP_A")
consequences <- missense_e
#minReadShift <- 10
#minMafDif <- 0.5

for (arg in cmd_args) {
  write(arg,stderr());
  arg <- strsplit(arg,"=")[[1]];
  if(arg[1] == "minMafDif") minMafDif <- as.numeric(arg[2]);
  if(arg[1] == "minReadShift") minReadShift <- as.numeric(arg[2]);
  if(arg[1] == "windows") windowsFile <- as.character(arg[2]);
  if(arg[1] == "vcf") vcfFileName <- as.character(arg[2]);
  if(arg[1] == "vep") vepFile <- as.character(arg[2]);
  if(arg[1] == "coding") consequences <- coding_e;
  if(arg[1] == "nonsense") consequences <- nonsense_e;
  
#  if(arg[1] == "missense") consequences <- missense_e;
}

tryCatch(
  {
    
        
    windows <- read.table(windowsFile)
    if(dim(windows)[[2]] == 1) {
      windows <- t(matrix(unlist(strsplit(as.character(windows[,1]),split="[:-]")),nrow=3))
      }
    windows <- as.data.frame(windows)
    colnames(windows) <- c("chr","start","end")
    win.gr <- GRanges(seqnames=windows$chr, ranges=IRanges(as.numeric(as.character(windows$start)), end=as.numeric(as.character(windows$end))))
    
    head <- scanVcfHeader(paste(vcfFileName,".gz", sep=""))
    write(paste(unlist(head)),stderr())

    vcftab <- TabixFile(paste(vcfFileName,".gz", sep=""), index = paste(vcfFileName, "gz.tbi", sep="."));
    vcftab

    myParam <- ScanVcfParam(which=win.gr, GENO=c("DP4","DP"),INFO=c("CSQ","FS"))
    wins.vcf.scan <- scanVcf(paste(vcfFileName,".gz", sep=""), param=myParam)

    #REPLACE NAs with 0,0,0,0
    wins.vcf.scan[[1]]$GENO$DP4[which(is.na(wins.vcf.scan[[1]]$GENO$DP4))] <- "0,0,0,0"
    
    write("unpacking",stderr())
    wins.vcf <- unpackVcf(wins.vcf.scan, vcftab, info=TRUE, geno=TRUE)
    
    wins.vcf <- getBafsMultiFile(wins.vcf)   
    wins.vcf <- getMaxBafContrast(wins.vcf)   
    wins.vcf <- getMaxReadShift(wins.vcf)

    if(exists("vepFile")) {
      genes <-  read.table(vepFile, comment.char="#")
      colnames(genes) <- c("Uploaded_variation","Location","Allele","Gene",
                           "Feature","Feature_type","Consequence","cDNA_position",
                           "CDS_position","Protein_position","Amino_acids","Codons",
                           "Existing_variation","Extra")
      var_ids <- as.data.frame(t(matrix(unlist(strsplit(as.character(genes$Uploaded_variation),split="[_]",)),nrow=3)))
      genes <- cbind(var_ids[,1:2], genes$Allele, genes$Gene, genes$Feature, genes$Consequence)
      colnames(genes) <- c("CHR","POS","ALLELE","GENE","FEATURE","CSQ")
      printValidVars(wins.vcf,minMC=minMafDif, minRS=minReadShift,CSQs=consequences, genes=genes)
    }
    else {    
      # printValidVars(wins.vcf,minMC=0.5,CSQs=missense_e)
      printValidVars(wins.vcf,minMC=minMafDif, minRS=minReadShift,CSQs=consequences)
    }
    #WRITE THE REST OF THIS SCRIPT!!!!!
    # print N if MBC >= 0.5... (or whatever) 
    # bafs <- wins.vcf[[i]]$GENO$BAF
    
    
    }, #end of try
    warnings = function(w) {w
                              print(w$message,stderr())},
    error = function(e) {
                              if (e$message == "no lines available in input"){
                                  warning(e$message)
                                  quit(status=0);
                                  }
                              e
                              },
         finally = { write("end",stderr()); }
  
)

