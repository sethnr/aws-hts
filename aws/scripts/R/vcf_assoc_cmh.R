# load required libraries
library("zoo")
library("ShortRead")
library("snpStats");
library("GGtools");
# library("VariantAnnotation")

########
# standard modules (shift to dedicated library)
########

writeVcfFromStdin <- function(filename) {
  vcfFile <- file(filename, "w")
  fin <- file("stdin")
  
  open(fin)
  while(length(line <- readLines(fin,n=1)) > 0) {
    #nb, scanVcf will flip if given trailing tabs, chomp everything. 
    write(chomp(line), vcfFile)
  }           
  close(fin)
  close(vcfFile)  
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
}

# read in table from STDIN
tryCatch(
  {
          vcfFileName <- "tmpvcf.vcf"
          writeVcfFromStdin(vcfFileName)  
          system(paste("./bgzip ",vcfFileName), wait=T)
          system(paste("./tabix -pvcf ",vcfFileName,".gz", sep=""), wait=T)
          system(paste("ls -l ",vcfFileName,"*  2>&1 ", sep=""))


  
#           vcfFileName <- "~/Scratch/theta/tmpvcf.vcf";
#	sessionInfo()
#	   testVcfFile <- file(vcfFileName)
#	testVcfFile
#	while(length(line <- readLines(testVcfFile)) > 0) {write(line,stderr())}
# window <- 10000; step <- 1000; maxDp<-20
           
	         vcf <- getVCF(paste(vcfFileName,".gz", sep=""))
           head <- scanVcfHeader(paste(vcfFileName,".gz", sep=""))
           write(paste(unlist(head)),stderr())
           
           
           
          vcf <- getBafsMultiFile(vcf)           
           
           # write top two lines & dimensions to STDERR
           # depth <- getDepthTable(depth)
           #if(fillBases) { depth <- fillInMissingBases(depth)}
           theta <- calcSWThetaPerCs(vcf, window=window, step=step, maxDp=max_pool_size)
           
           write.table(theta,stdout(), row.names=F, col.names=F)
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

