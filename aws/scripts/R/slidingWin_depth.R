# load required libraries
library("zoo")

getDepthTable <- function(depth, sort) {
  write("in:",stderr())
  write.table(depth[1:4,],stderr(), row.names=F, col.names=T)
  #direct output from aggregate -> make into 3-col output
  if(dim(depth)[[2]] == 2) {
#    write("cols = 2",stderr())
#    write.table(depth[1:2,],stderr(), row.names=F, col.names=T)
    pos_chr <- strsplit(sapply(depth, function(x) {gsub("LongValueSum:","",x)}),"\\.")
    chrs <- unlist(pos_chr)[1:dim(depth)[[1]]*2-1]
    chrs <- as.numeric(chrs)
    poss <- unlist(pos_chr)[1:dim(depth)[[1]]*2]
    poss <- as.numeric(poss)  
    
    depth <- as.data.frame(cbind(chr=chrs,pos=poss,dep=depth[,2]))
#    names(depth) <- c("chr","pos","dep")
#    depth <- depth[with(depth, order(chr, pos)), ]
  }
  #split aggregated output (chr,pos,dep)
  if(dim(depth)[[2]] == 3) {
#    write("cols = 3",stderr())
#    write.table(depth[1:2,],stderr(), row.names=F, col.names=T)
    names(depth) <- c("chr","pos","dep")
    depth <- depth[with(depth, order(chr, pos)), ]
  }
  #granges output (i.e. output of this script)
  if(dim(depth)[[2]] == 10) {
#    write("cols = 10",stderr())
#    write.table(depth[1:2,],stderr(), row.names=F, col.names=T)
    names(depth) <- c("chr","start","end","length","strand","name","num.mark","dep","startRow","endRow")
    depth <- depth[with(depth, order(chr, start)), ]
  }

  write("out:",stderr())
  write.table(depth[1:4,],stderr(), row.names=F, col.names=T)
depth
}

fillInMissingBases <- function(depth, step) {
  for (chr in unique(depth$chr)) {chrDep <- depth[which(depth$chr == chr),];
                                  depth <- depth[which(depth$chr != chr),];
                                  fl <- floor((min(chrDep$pos)/step))*step;
                                  ce <- ceiling((max(chrDep$pos)/step))*step;
                                  posns <- fl:ce
                                  newDep <- data.frame(rep(chr,length(posns)),posns);
                                  names(newDep) <- c("chr","pos");
                                  newDep <- merge(newDep,chrDep,by=c("chr","pos"),all=T);
                                  depth <- rbind(depth,newDep);
                                }
  depth
}

calcSWPerCs <- function(depth, window, step) {
  for (chr in unique(depth$chr)) {
    chrDep <- depth[which(depth$chr == chr),];

    write(paste("processing chunk ",chr," by ",step),stderr())
    write.table(dim(chrDep),stderr(), row.names=T, col.names=T)
    write.table(summary(chrDep),stderr(), row.names=F, col.names=F)


    depthW <- rollapply(chrDep$dep, width=window, by=step, na.rm=T, FUN=mean)
    startW <- rollapply(chrDep$pos, width=window, by=step, na.rm=T, FUN=min)
    endW <- rollapply(chrDep$pos, width=window, by=step, na.rm=T, FUN=max)
    chrDepSW <- cbind(rep(chr,length(endW)),startW,endW,depthW);
    names(chrDepSW) <- c("chr","start","end","depth")
    if(exists("depthSW")) {depthSW <- rbind(depthSW,chrDepSW);}
    else {depthSW <- chrDepSW;
          names(depthSW) <- c("chr","start","end","depth")}
    
  }
  depthSW
}


# parse arguments from command line
minseg <- 0
cmd_args = commandArgs(trailingOnly = T);
for (arg in cmd_args) {
  write(arg,stderr());
  arg <- strsplit(arg,"=")[[1]];
  if(arg[1] == "window") window <- as.numeric(arg[2]);
  if(arg[1] == "step") step <- as.numeric(arg[2]);
  if(arg[1] == "fillBases") fillBases <- arg[2];
}
if(!exists("fillBases")) {fillBases <- FALSE}

# read in table from STDIN
tryCatch(
         {
           depth <- read.table("stdin", header=F);
           # write top two lines & dimensions to STDERR
           depth <- getDepthTable(depth)
           #if(fillBases) { depth <- fillInMissingBases(depth)}
           #depth <- calcSWPerCs(depth, window=window, step=step)
           
           #write.table(depth,stdout(), row.names=F, col.names=F)
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

