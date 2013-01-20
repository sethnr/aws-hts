getWindows <- function(ld, window, step) {
#  outcols <-colnames(ld)
#  outcols[2:3] <- c("start","end")
  windows <- data.frame()
  for (chr in unique(ld$chr)) {
   ldc <- ld[which(ld$chr == chr),]
   min <- floor((min(as.numeric(as.vector(ldc$SNP1)))/step))*step;  
   max <- ceiling((max(as.numeric(as.vector(ldc$SNP2)))/step))*step; 
   w_start <- seq(min,(max-window),by=step)
   w_end <- seq((min+window),max,by=step)
   windows <- rbind(windows,data.frame("chr"=rep(chr,length(w_start)),w_start,w_end));
  }
  windows
}

windowMeans <- function(x,m) {
  chr <- x[1]
  w_start <- as.numeric(x[2])
  w_end <- as.numeric(x[3])
  
  cols <- dim(m)[[2]]
  rows <- 0
  #figure this shit out!!!
  wm <- m[which(m$chr==chr),]
  wm <- wm[which(wm[2] >= w_start),]
  wm <- wm[which(wm[3] <= w_end),]
  rows <- dim(wm)[[1]]
  csums <- colSums(wm[,4:cols])
  outline <- c(chr,"START"=w_start, "END"=w_end, csums / rows)
  write(cat(outline,sep="\t"), stdout())
  }



# parse arguments from command line
step <- 2000
window <- 10000
cmd_args = commandArgs(trailingOnly = T);
for (arg in cmd_args) {
  write(arg,stderr());
  arg <- strsplit(arg,"=")[[1]];
  if(arg[1] == "window") window <- as.numeric(arg[2]);
  if(arg[1] == "step") step <- as.numeric(arg[2]);
}

#setwd("/Users/seth/Scratch/deepseq_LD")
#file <- "Fd03.ld"

# read in table from STDIN
ld <- read.table(file("stdin"))
colnames(ld) <- c("chr","SNP1", "SNP2",
                  #Number of pairs observed with: 
                  "x_11", "x_12", "x_21", "x_22",
                  #Estimate for allele frequency of allele: 
                  "fA","fB",
                  #read depth SNP:
                  "d1","d2",
                  #Intersecting read depth
                  "dI",
                  #MLE estimate, low, actual, high:                  
                  "MLE_low", "R2_MLE","MLE_high",
                  #Direction Computation R2
                  "R2_dir",
                  #alleles: 
                  "A","B","a","b")

ld <- ld[1:(dim(ld)[[2]]-4)]
windows <- getWindows(ld, window, step)


outcols <-colnames(ld)
outcols[2:3] <- c("start","end")
write(cat(outcols,sep="\t"), stdout())
apply(windows,1,windowMeans,m=ld)


# SEGS DO NOT FIND SEPARATE REGIONS
#ld_gr <- GRanges(seqnames=Rle(ld$chr), ranges=IRanges(ld$SNP1, end=ld$SNP2), depth=ld$R2_MLE)
#ld_segs <- fastseg(ld_gr)
#ld_segs <- as.data.frame(ld_segs)[,c("seqnames","start","end","seg.mean")]
#write.table(ld_segs,stderr())