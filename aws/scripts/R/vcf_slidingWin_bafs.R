# load required libraries
#library("zoo")
library("ShortRead")
library("snpStats");
library("GGtools");
library("VariantAnnotation")

########
# standard modules (shift to dedicated library)
########

getWindows <- function(window, step, chr, min, max) {
  w_start <- seq(min,(max-window),by=step)
  w_end <- seq((min+window),max,by=step)
  windows <- data.frame("chr"=rep(chr,length(w_start)),w_start,w_end);
#  windows[1,]
  windows
  }

calcSWThetaPerCs <- function(vcf, window, step, maxDp, chunk) {
  baf <- vcf[[1]]$GENO$BAF
  baf[which(is.na(baf))] <- 0 
  baf_dp <- as.data.frame(cbind(vcf[[1]]$CHR,vcf[[1]]$POS,baf,vcf[[1]]$GENO$DP))
  # colnames(baf_dp)[1:2] <- c("chr","pos")
  colnames(baf_dp) <- c("chr","pos",sub(".vcf",".BAF",colnames(vcf[[1]]$GENO$GT)),sub(".vcf",".DP",colnames(vcf[[1]]$GENO$GT)))
  
  
#  if(exists("maxDp")) {baf_dp <- limitDps(baf_dp,maxDp)}
  
  for (chr in unique(baf_dp$chr)) {
    chrBafDp <- baf_dp[which(baf_dp$chr == chr),];

#    write.table(summary(chrBafDp),stderr(), row.names=F, col.names=F)
    
    write("calc theta ",stderr())
    write.table(head(chrBafDp, n=1L),col.names=F,row.names=F,stderr())
    write.table(tail(chrBafDp, n=1L),col.names=F,row.names=F,stderr())
    fl <- floor((min(as.numeric(as.vector(chrBafDp$pos)))/step))*step;
    ce <- ceiling((max(as.numeric(as.vector(chrBafDp$pos)))/step))*step;
    windows <- getWindows(window, step, chr, fl, ce)
    
#    rollapply(chrBafDp, width=window, by=step, by.column=F, FUN=calcThetaBAFs)
    chrBafDp <- as.data.frame(chrBafDp)
#    apply(windows,1,runFuncOnMatrixWin,m=chrBafDp, func=calcThetaBAFs)
    apply(windows,1,runFuncOnMatrixWin,m=chrBafDp, func=calcHpBAFs)
    write("rollApply finished",stderr())
    }
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
  colnames(bafs) <- colnames(vcf[[1]]$GENO$DP4)
  vcf[[1]]$GENO$BAF <- bafs;
  vcf
}

permuteVCF <- function(vcf, n) {
  noSamples <- dim(vcf[[1]]$GENO$DP)[2];
  noVars <- dim(vcf[[1]]$GENO$DP)[1];
  noVals <- length(vcf[[1]]$GENO$DP);
  
  index <- matrix(1:noVals,ncol=noSamples)
  p_index <- t(apply(index,1,FUN=sample)) 
  newcolnames <- paste(colnames(vcf[[1]]$GENO$DP),".",n,sep="")
  vcf[[1]]$GENO$DP <- matrix(vcf[[1]]$GENO$DP[p_index],ncol=noSamples) 
  colnames(vcf[[1]]$GENO$DP) <- newcolnames
  vcf[[1]]$GENO$BAF <- matrix(vcf[[1]]$GENO$BAF[p_index],ncol=noSamples) 
  colnames(vcf[[1]]$GENO$BAF) <- newcolnames
  
  vcf
}



runFuncOnMatrixWin <- function(x,m,func) {
  chr <- x[1]
  
  #strip nas from matrix
  #WHY THE HELL ARE THEY HERE???
  
  w_start <- as.numeric(x[2])
  w_end <- as.numeric(x[3])
#  write(paste(" ** ",chr,":",w_start,"-",w_end," ** ", sep=""),stderr())
  posns <- as.numeric(as.vector(m[[2]]))
  m <- m[which(m[[1]]==chr),]
  m <- m[which(posns >= w_start),]
  m <- m[which(posns <= w_end),]
  #m <- m[which(!is.na(posns)),]
  m <- m[which(!is.na(m[,2])),]
  m <- as.matrix(m)
  
#  write.table(m,stderr())
  func(m, w_start, w_end)
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


calcHp <- function(dp,maf) {
  # if given BAF, silently convert to MAF (flip any majors top minor)
  maf[which(maf > 0.5)] <- 1-maf[which(maf > 0.5)]
  # major & non-major allele counts
  majDp <- (1-maf) * dp
  minDp <- maf * dp;
  
  Hp <- (2 * sum(majDp) * sum(minDp)) / (sum(dp)^2)
  Hp
}    


calcHpBAFs <- function(x, w_start, w_end) {
  # x = BAF/DP matrix
  if(!is.matrix(x)) {write("I was expecting a matrix", stderr());
  }
  
  noSamples <- (dim(x)[2]-2) /2;
  Hps <- vector()
  posns <- as.numeric(x[,2])

  #prep output line
  if(exists("w_start") && exists("w_end")) {
    Hps <- c(Hps,x[1,1],w_start,w_end)  }
  else {
    Hps <- c(Hps,x[1,1],min(posns), max(posns))  }

  sumDp <- rep.int(0,length(posns))
  bDp <- rep.int(0,length(posns))
  for(s in 1:noSamples) {
    dp <- as.numeric(x[,s+2+noSamples])
    #nb can set NAs in depth to 0, as we're working with totals will make no difference 
    #(unless all are NA in which case it'll probably flip the fuck out)
    dp[which(is.na(dp))] <- 0 

    baf <- as.numeric(x[,s+2])
    baf[which(is.na(baf))] <- 0
    
    # save sums for calculation of total
    bDp <- (bDp + (baf * dp))
    sumDp <- sumDp + dp
    
    Hp <- calcHp(dp,baf)
    Hps <- c(Hps,Hp)
    rm(Hp,dp,baf)
    gc()
  }  
  sumBaf <- bDp / sumDp
  sumHp <- calcHp(sumDp,sumBaf)
  Hps <- c(Hps,sumHp)
  
  rm(sumBaf, sumHp, sumDp)  

  write(cat(Hps,sep="\t"),stdout())
  rm(Hps)
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
input_type <- "vcf"
permute <- -1 #null value (no permutation)
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
  if(arg[1] == "permute") permute <- as.numeric(arg[2]);
  if(arg[1] == "input") input_type <- arg[2];
  if(arg[1] == "input_type") input_type <- arg[2];
}
if(!exists("max_pool_size")) {max_pool_size<-100}

write("all configs read and that",stderr())

write("tar zxvf *.tar.gz",stderr())
system("ls -l 1>&2",wait=T)
system("tar zxvf *.tar.gz",wait=T);  
system("tar zxvf *.tgz",wait=T);  
# system("ls -l 1>&2",wait=T) 


tryCatch(
  {

vcfFileNames <- vector()
if(input_type == "vcf") {
  write("reading VCF from STDIN",stderr())
  vcfFileName <- "tmpvcf.vcf"
  # read in table from STDIN
  fin <- file("stdin")
  open(fin)
    
  vcfFile <- file(vcfFileName, "w")
  while(length(line <- readLines(fin,n=1)) > 0) {
  #nb, scanVcf will flip if given trailing tabs, chomp everything. 
    write(chomp(line), vcfFile)
    }           
  close(fin)
  close(vcfFile)
  vcfFileNames <- c(vcfFileNames, vcfFileName);
  vcfIndex <- cbind(vcfFileNames,rep(permute,length(vcfFileNames)))
}


if (input_type == "index") {
  vcfIndex <- read.table(file("stdin"))
  noFiles <- (length(vcfIndex)/2)
  write(paste("fetching ",noFiles,"files"),stderr())

  vcfLocalFiles <- vector()
  vcfPerms <- vcfIndex[,2]
  for (i in 1:noFiles) {
    vcfFileName <- vcfIndex[i,1]
    write(paste(vcfFileName," is.character = ",is.character(vcfFileName)),stderr())
    write(paste(vcfFileName," is.vector = ",is.vector(vcfFileName)),stderr()) 
    vcfBaseName <-  basename(as.character(vcfFileName))
    if ((substr(vcfFileName,1,5) == "s3://") || (substr(vcfFileName,1,6) == "s3n://")) {     
      write(paste("s3cmd  --config ./s3cmd.config get ",vcfFileName," ",vcfBaseName),stderr())
      system(paste("s3cmd  --config ./s3cmd.config get ",vcfFileName," ",vcfBaseName),wait=T);  
      }
    else {
      write(paste("wget ",vcfFileName,"./"),stderr())
      system(paste("wget ",vcfFileName,"./"),wait=T);  
      }
    vcfLocalFiles <- c(vcfLocalFiles,as.character(vcfBaseName))
    
    write(paste(vcfIndex[i,1],vcfBaseName),stderr())
    }
vcfIndex <- cbind(vcfLocalFiles,vcfPerms)
}
write("files fetched from S3/elsewhere",stderr())



for (i in 1:(length(vcfIndex)/2)) {
  write(paste("processing file",vcfIndex[i,1],"perms",vcfIndex[i,2]),stderr())
  vcfFileName <- vcfIndex[i,1]
  write(paste("zipping and indexing ",vcfFileName),stderr())
  system(paste("./bgzip ",vcfFileName), wait=T)
  system(paste("./tabix -pvcf ",vcfFileName,".gz", sep=""), wait=T)
}
system("ls -l * 1>&2", wait=T)

#sessionInfo()

#tryCatch(   # MOVED TRYCATCH ABOVE TO CATCH EMPTY-INDEX-FILES (make separate tomorrow)
#  {           
  #           vcfFileName <- "~/Scratch/theta/3R.6492534.31999987.vcf";
  #           window <- 10000; step <- 1000; maxDp<-20; chunk <- 500000
  
        for (i in 1:(length(vcfIndex)/2)) {
        vcfFileName <- vcfIndex[i,1]
        permute <- vcfIndex[i,2]}
  
    
    	  vcf <- getVCF(paste(vcfFileName,".gz", sep=""))
        head <- scanVcfHeader(paste(vcfFileName,".gz", sep=""))
        write(paste(unlist(head)),stderr())
    
        vcf <- getBafsMultiFile(vcf)

        baf <- vcf[[1]]$GENO$BAF
        dp <- vcf[[1]]$GENO$DP
        write.table(summary(vcf[[1]]),stderr())
        totalBAF <- apply((baf * dp),1,FUN=sum) / apply(dp,1,FUN=sum)
        baf <- cbind(vcf[[1]]$CHROM,vcf[[1]]$POS,baf,totalBAF)
        write.table(baf,col.names=F, row.names=F, quote=F, sep="\t",stdout())
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

