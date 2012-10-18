# load required libraries
library("zoo")

getDepthTable <- function(depth, sort) {
  #direct output from aggregate -> make into 3-col output
  if(dim(depth)[[2]] == 2) {
    pos_feat <- strsplit(sapply(depth, function(x) {gsub("LongValueSum:","",x)}),"\\.")
    feats <- unlist(pos_feat)[1:dim(depth)[[1]]*2-1]
    feats <- as.numeric(feats)
    poss <- unlist(pos_feat)[1:dim(depth)[[1]]*2]
    poss <- as.numeric(poss)  
    
    depth <- as.data.frame(cbind(feat=feats,pos=poss,dep=depth[,2]))
    names(depth) <- c("feat","pos","dep")
    depth <- depth[with(depth, order(feat, pos)), ]
  }
  #split aggregated output (feat,pos,dep)
  if(dim(depth)[[2]] == 3) {
    names(depth) <- c("feat","pos","dep")
    depth <- depth[with(depth, order(feat, pos)), ]
  }
  #granges output (i.e. output of this script)
  if(dim(depth)[[2]] == 10) {
    names(depth) <- c("feat","start","end","length","strand","name","num.mark","dep","startRow","endRow")
    depth <- depth[with(depth, order(feat, start)), ]
  }  
depth
}

summaryPerFeat <- function(depth) {
  for (feat in unique(depth$feat)) {
    featDep <- depth[which(depth$feat == feat),];   
#    write.table(dim(featDep),stderr(), row.names=T, col.names=T)
    if(exists("featSum")) {featSum <- rbind(featSum,c(feat,mean(featDep$dep),sd(featDep$dep)));}
    else {featSum <- c(feat,mean(featDep$dep),sd(featDep$dep));
          names(featSum) <- c("feat","mean","depth")}
    
  }
  featSum
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
           write.table(depth[1:2,],stderr(), row.names=F, col.names=T)
           write(dim(depth),stderr())
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
             #if(fillBases) { depth <- fillInMissingBases(depth,step)}
             #depth <- calcSWPerCs(depth, window=window, step=step)
             #write.table(depth,stdout(), row.names=F, col.names=F)
             featSum <- summaryPerFeat(depth);
             write.table(featSum,stdout(), row.names=F, col.names=F, quote=F)
           }
         }
)

