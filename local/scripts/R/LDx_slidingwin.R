library("GGtools");

library("zoo");

window <- 1000
step <- 200
cmd_args = commandArgs(trailingOnly = T);
for (arg in cmd_args) {
  write(arg,stderr());
  arg <- strsplit(arg,"=")[[1]];
  if(arg[1] == "file") ld_pairs_file <- as.character(arg[2]);
  if(arg[1] == "window") window <- as.numeric(arg[2]);
  if(arg[1] == "step") step <- as.numeric(arg[2]);
}
LD_pairs <- read.table(ld_pairs_file, header=T);

noCols <- length(colnames(LD_pairs))

if(colnames(LD_pairs)[[2]] == "pos") {
  write(paste(c("chr","start","end",colnames(LD_pairs)[3:noCols]),sep="\t",collapse="\t"),stdout())
  } else {
  write(paste(colnames(LD_pairs),sep="\t",collapse="\t"),stdout())
}


for (chr in unique(LD_pairs$chr)) {
  LDc_pairs <- LD_pairs[which(LD_pairs$chr == chr),]
  write(paste(chr,dim(LDc_pairs)[[1]]),stderr())
  if (dim(LDc_pairs)[[1]] > window) {
    if(colnames(LD_pairs)[[2]] == "pos") {
      chrs <- rollapply(as.vector(LDc_pairs$chr),window,FUN=min,by=step)
      starts <- rollapply(LDc_pairs$pos,window,FUN=min,by=step)
      ends <- rollapply(LDc_pairs$pos,window,FUN=max,by=step)
      posTable <- cbind(chrs,starts, ends)
      index<-3
    }
    else {
      chrs <- rollapply(as.vector(LDc_pairs$chr),window,FUN=min,by=step)
      starts <- rollapply(LDc_pairs$start,window,FUN=min,by=step)
      ends <- rollapply(LDc_pairs$end,window,FUN=max,by=step)
      index<-4
      posTable <- cbind(chrs,starts,ends)    
    }
  
    for (i in index:length(colnames(LD_pairs))) {
      vals <- rollapply(LDc_pairs[,i],window,FUN=mean,by=step)
      valName = colnames(LD_pairs)[[i]]
      posTable <- cbind(posTable,valName=vals)
      }
    rm(LDc_pairs);
    write.table(posTable,quote=F, col.names=F, row.names=F, stdout());  
    }
  }