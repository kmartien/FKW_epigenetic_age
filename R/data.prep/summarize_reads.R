library(readr)
library(seqinr)

summarize.reads <- function(min.cov) {

  description <- paste0("min.cov.",min.cov)
  
  bs.refs <- read.fasta("data_raw/Pcra.bisulfite.amplicon.wPrimers.fasta")
  orig.refs <- read.fasta("data_raw/Pcra.unconverted.amplicon.wPrimer.fasta")
  primers <- read.fasta("data_raw/Pcra.epi.primers.fasta")
  amplicons <- list("TET2","FLJ0945","DIRAS3","VGF","HMG20B","KCNC4","PRDM12","GRIA2")
  
  #determine the locations within the references that aren't part of the primers
  target.positions <- sapply(names(bs.refs), function(p){
    target.start <- length(primers[[which(names(primers)==paste(p,"_F1", sep=""))]]) +1
    target.end <- length(bs.refs[[which(names(bs.refs)==p)]]) - length(primers[[which(names(primers)==paste(p,"_R1", sep=""))]])
    return(c(start=target.start, end=target.end))
  })
  
  pileup.path <- "data_raw/MS15.galaxy.alignments/Aligned.to.targets.w.primers/pileups"
  fnames <- list.files(path = pileup.path, 
                       pattern="\\.fastq.tabular.csv$")
  num.samps <- length(fnames)
  
  plp.all <- do.call("rbind",lapply(fnames, function(f){
    read.csv(paste0(pileup.path,"/",f)) %>% as.data.frame
  }))
  
  id <- do.call("c",lapply(plp.all$fname, function(f){
    strsplit(f,"_",fixed=TRUE)[[1]][1]
  }))
  plp.all$fname <- id
  names(plp.all)[c(1,3)] <- c("id","position")

  meth.dat <- lapply(amplicons, function(amp){
    sub.plp <- subset(plp.all,chrom==amp)
    bs.ref.seq <- bs.refs[[which(names(bs.refs)==amp)]]
    orig.ref.seq <- orig.refs[[which(names(orig.refs)==amp)]]
    target.pos <- target.positions[,which(colnames(target.positions)==amp)]
    calc_methylation(sub.plp,bs.ref.seq, orig.ref.seq, target.pos, min.cov)
  })
  names(meth.dat) <- amplicons
  
  save(meth.dat, file = paste0("data/meth.dat.",description,".Rdata"))
  return(meth.dat)
}
