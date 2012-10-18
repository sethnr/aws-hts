# load required libraries
library("fastseg")

getGRanges <- function(depth) {
  #direct output from aggregate -> make into 3-col output
  if(dim(depth)[[2]] == 2) {
    pos_chr <- strsplit(sapply(depth, function(x) {gsub("LongValueSum:","",x)}),"\\.")
    chrs <- unlist(pos_chr)[1:dim(depth)[[1]]*2-1]
    chrs <- as.numeric(chrs)
    poss <- unlist(pos_chr)[1:dim(depth)[[1]]*2]
    poss <- as.numeric(poss)  
    
    depth <- as.data.frame(cbind(chr=chrs,pos=poss,dep=depth[,2]))
    names(depth) <- c("chr","pos","dep")
    depth <- depth[with(depth, order(chr, pos)), ]
  }
  #split aggregated output (chr,pos,dep)
  if(dim(depth)[[2]] == 3) {
    names(depth) <- c("chr","pos","dep")
    depth <- depth[with(depth, order(chr, pos)), ]
    gr <- GRanges(seqnames=depth$chr,
                  ranges=IRanges(depth$pos, end=depth$pos),depth=as.double(depth$dep))
  }
  #granges output (i.e. output of this script)
  if(dim(depth)[[2]] == 10) {
    names(depth) <- c("chr","start","end","length","strand","name","num.mark","dep","startRow","endRow")
    depth <- depth[with(depth, order(chr, start)), ]
    gr <- GRanges(seqnames=depth$chr,
                  ranges=IRanges(depth$start, end=depth$end),depth=as.double(depth$dep))
  }
gr
}




# parse arguments from command line
minseg <- 0
cmd_args = commandArgs(trailingOnly = T);
for (arg in cmd_args) {
  arg <- strsplit(arg,"=")[[1]];
  write(paste(arg[1]," -> ",arg[2]),stderr());
  if(arg[1] == "minSeg") minseg <- arg[2];
  }

# read in table from STDIN
tryCatch(
         {
         depth <- read.table("stdin", header=F);
         # write top two lines & dimensions to STDERR
         write.table(depth[1:2,],stderr(), row.names=F, col.names=T)
         write(dim(depth),stderr())

         gr <- getGRanges(depth)
         segs <- fastseg(gr,minSeg=minseg)
         segs <- as.data.frame(segs)
         write.table(segs, stdout(), row.names=F, col.names=F)
         
         },
         warnings = function(w) {w
                                print(w$message,stderr())
                               },
         error = function(e) {
           if (e$message == "no lines available in input"){
            warning(e$message)
            quit(status=0);
           }
           e}
         ,
         finally = {}
)

