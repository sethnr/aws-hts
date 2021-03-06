# load required libraries
library("zoo")
library("ShortRead")
library("snpStats");
library("GGtools");
# library("VariantAnnotation")

########
# standard modules (shift to dedicated library)
########

calcSWThetaPerCs <- function(vcf, window, step, maxDp, chunk) {
  baf <- vcf[[1]]$GENO$BAF
  baf[which(is.na(baf))] <- 0 
  baf_dp <- as.data.frame(cbind(vcf[[1]]$CHR,vcf[[1]]$POS,baf,vcf[[1]]$GENO$DP))
  # colnames(baf_dp)[1:2] <- c("chr","pos")
  colnames(baf_dp) <- c("chr","pos",sub(".vcf",".BAF",colnames(vcf[[1]]$GENO$GT)),sub(".vcf",".DP",colnames(vcf[[1]]$GENO$GT)))
  
  
#  if(exists("maxDp")) {baf_dp <- limitDps(baf_dp,maxDp)}
  baf_dp <- fillInMissingBases(baf_dp, step, chunk)
  
  for (chr in unique(baf_dp$chr)) {
    chrBafDp <- baf_dp[which(baf_dp$chr == chr),];

#    write.table(summary(chrBafDp),stderr(), row.names=F, col.names=F)
    
    write("calc theta ",stderr())
    write.table(head(chrBafDp, n=1L),col.names=F,row.names=F,stderr())
    write.table(tail(chrBafDp, n=1L),col.names=F,row.names=F,stderr())

    
    rollapply(chrBafDp, width=window, by=step, by.column=F, FUN=calcThetaBAFs)
    write("rollApply finished",stderr())
#    names(chrDepSW)[1:3] <- c("chr","start","end")
#    if(exists("bafDpSW")) {bafDpSW <- rbind(bafDpSW,chrDepSW);}
#    else {bafDpSW <- chrDepSW;
#          names(bafDpSW)[1:3] <- c("chr","start","end")}    
  }
#  bafDpSW
}

fillInMissingBasesMC <- function(depth, step, chunk) {
#  returnDep <- Matrix(ncol=ncol(depth), sparse=T)
  for (chr in unique(depth$chr)) {chrDep <- depth[which(depth$chr == chr),];
                                  depth <- depth[which(depth$chr != chr),];
                                  
#                                  fl <- floor((min(as.numeric(as.vector(chrDep$pos)))/step))*step;
#                                  last <- fl;
#                                  posns <- vector();
#                                  for (pos in as.numeric(as.vector(chrDep$pos))) {
#                                    if (pos > last+chunk) {
#                                      ce <- ceiling((last/step))*step;  
#                                      posns <- c(posns,fl:ce);
#                                      write(paste(fl,ce,length(posns)),stderr())
#                                      fl <- floor(pos/step)*step;
#                                      }
#                                    last <- pos
#                                    }
#                                  ce <- ceiling((last/step))*step;
#                                  posns <- c(posns,fl:ce);
#                                  write(paste(fl,ce,length(posns)),stderr())
                                  
                                  fl <- floor((min(as.numeric(as.vector(chrDep$pos)))/step))*step;
                                  ce <- ceiling((max(as.numeric(as.vector(chrDep$pos)))/step))*step;
                                  
                                  posns <- fl:ce
                                  newDep <- data.frame(rep(chr,length(posns)),posns);
                                  names(newDep) <- c("chr","pos");
                                  newDep <- merge(newDep,chrDep,by=c("chr","pos"),all=T);
                                  
#                                  returnDep <- rbind(returnDep,newDep);
                                  
                                  if(exists("returnDep")) {returnDep <- rbind(returnDep,newDep);}
                                  else {returnDep <- newDep;}
  }
  returnDep
}

fillInMissingBases <- function(depth, step, chunk) {
  fl <- floor((min(as.numeric(as.vector(depth$pos)))/step))*step;
  ce <- ceiling((max(as.numeric(as.vector(depth$pos)))/step))*step;                                  
  posns <- fl:ce
  if (length(unique(depth$chr)) > 1) {# throw error 
    stop(paste("multiple chrs found",unique(depth$chr),":",fl,"-",ce))}
  chr <- unique(depth$chr)[[1]]
  newDep <- data.frame(rep(chr,length(posns)),posns);
  names(newDep) <- c("chr","pos");
  newDep <- merge(newDep,depth,by=c("chr","pos"),all=T);
  newDep      
}


limitDps <- function(x, maxDp) {
  noSamples <- (dim(x)[2]-2) /2;
  for (i in 1:noSamples) {
    j <- i+2+noSamples
    dp <- as.numeric(x[,j])
    dp[which(dp > maxDp)] <- maxDp
    x[,j] <- dp
    }
  x
  }



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

calcThetaBAFs <- function(x) {
  # x = BAF/DP matrix
  noSamples <- (dim(x)[2]-2) /2;
  thetas <- vector()
#  thetas <- c(thetas, )
  posns <- as.numeric(x[,2])
#  thetas <- c(thetas,x[1,1],min(x[,2]), max(x[,2]))  
  thetas <- c(thetas,x[1,1],min(posns), max(posns))  
  
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

nthHarmonic <- function(n) { 
    h<-0; 
    if(!is.na(n)) {for (i in 1:n)
        { h<-(h+(1/i));}
                  } 
    h; 
    }


chomp <- function (x) sub("\\s+$", "", x)







# parse arguments from command line
cmd_args = commandArgs(trailingOnly = T);
for (arg in cmd_args) {
  write(arg,stderr());
  arg <- strsplit(arg,"=")[[1]];
  if(arg[1] == "window") window <- as.numeric(arg[2]);
  if(arg[1] == "step") step <- as.numeric(arg[2]);
  if(arg[1] == "pool_size") max_pool_size <- as.numeric(arg[2]);
  if(arg[1] == "max_pool_size") max_pool_size <- as.numeric(arg[2]);
  if(arg[1] == "max_size") max_pool_size <- as.numeric(arg[2]);
  if(arg[1] == "chunk") chunk <- as.numeric(arg[2]);
  if(arg[1] == "chunk_size") chunk <- as.numeric(arg[2]);
}
if(!exists("max_pool_size")) {max_pool_size<-100}

# read in table from STDIN
           vcfFileName <- "tmpvcf.vcf"
           vcfFile <- file(vcfFileName, "w")
           fin <- file("stdin")

           open(fin)
           while(length(line <- readLines(fin,n=1)) > 0) {
            #nb, scanVcf will flip if given trailing tabs, chomp everything. 
             write(chomp(line), vcfFile)
             }           
           close(fin)
           close(vcfFile)

#sessionInfo()


tryCatch(
{
  
           
           system(paste("./bgzip ",vcfFileName), wait=T)
           system(paste("./tabix -pvcf ",vcfFileName,".gz", sep=""), wait=T)
           system(paste("ls -l ",vcfFileName,"*", sep=""))


  
#           vcfFileName <- "~/Scratch/theta/tmpvcf.vcf";
#	sessionInfo()
#	   testVcfFile <- file(vcfFileName)
#	testVcfFile
#	while(length(line <- readLines(testVcfFile)) > 0) {write(line,stderr())}
#           vcfFileName <- "~/Scratch/theta/3R.6492534.31999987.vcf";
#           window <- 10000; step <- 1000; maxDp<-20; chunk <- 500000
           
	         vcf <- getVCF(paste(vcfFileName,".gz", sep=""))
           head <- scanVcfHeader(paste(vcfFileName,".gz", sep=""))
           write(paste(unlist(head)),stderr())
           
           
           
          vcf <- getBafsMultiFile(vcf)           
           
           # write top two lines & dimensions to STDERR
           # depth <- getDepthTable(depth)
           #if(fillBases) { depth <- fillInMissingBases(depth)}
           # theta <- calcSWThetaPerCs(vcf, window=window, step=step, maxDp=max_pool_size)
           calcSWThetaPerCs(vcf, window=window, step=step, maxDp=max_pool_size, chunk=chunk)
           
#           write.table(theta,stdout(), row.names=F, col.names=F)
         },
         warnings = function(w) {w
                                 print(w$message,stderr())},
         error = function(e) {
           if (e$message == "no lines available in input"){
             warning(e$message)
             quit(status=0);
           }
           e
         },
         finally = {
           if(exists("depth")) {
             #nb depth is calc'd per chunk (or CS if cs is in first col)
             if(fillBases) { depth <- fillInMissingBases(depth,step)}
             depth <- calcSWPerCs(depth, window=window, step=step)
             write.table(depth,stdout(), row.names=F, col.names=F)
           }
         }
)

