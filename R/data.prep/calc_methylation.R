calc_methylation <- function(plp, bs.ref, orig.ref, target.pos, min.cov = 100) {
  
  stopifnot(require(dplyr))
  stopifnot(require(tidyr))
  stopifnot(require(readr))
  stopifnot(require(seqinr))
  
  CpGs <- which(bs.ref=="c")
  all.Cs <- which(orig.ref=="c")
  
  #eliminate positions that are in the primer sequence
  CpGs <- CpGs[which(CpGs > target.pos[1] & CpGs < target.pos[2])]
  all.Cs <- all.Cs[which(all.Cs > target.pos[1] & all.Cs < target.pos[2])]
  
  base.cols <- c("A", "C", "G", "T", "N", "deletion","insertion")
  max.pos <- max(plp$position)
  
  # extract base frequencies from proportions
  base.freq <- round(plp[, base.cols] * plp$cov, 0)
  base.freq$position <- plp$position
  base.freq$id <- plp$id
  base.freq$coverage <- plp$cov
  
  # calculate % methylation (=Cs) at CpG sites and conversion efficiency at non-CpG Cs
  C.base.freq <- subset(base.freq,position %in% all.Cs)
  C.base.freq$non.CT <- rowSums(C.base.freq[,c(1,3,5,6)])
  C.summary <- lapply(all.Cs,function(x){
    coverage <- data.frame(subset(C.base.freq, position==x, select=c(id,coverage)))
    freq.meth <- data.frame(subset(C.base.freq, position==x, select=c(id,C)))
    freq.meth[which(coverage[,2] < min.cov),2] <- NA
    errors <- data.frame(subset(C.base.freq, position==x, select=c(id,non.CT)))
    names(coverage) <- names(freq.meth) <-names(errors) <- c("id", x)
    list(coverage=coverage, freq.meth=freq.meth, errors=errors)
  })
  names(C.summary) <- all.Cs
  
  CpG.summary <- C.summary[which(names(C.summary) %in% CpGs)]
  other.Cs.summary <- C.summary[-which(names(C.summary) %in% CpGs)]
  
  CpG.sum.df.list <- lapply(1:3, function(var){
    tab <- Reduce(function(x,y) full_join(x,y), lapply(CpG.summary, function(x){x[[var]]}))
    arrange(tab,id)
  })
  names(CpG.sum.df.list) <- c("coverage","freq.meth","errors")
  CpG.sum.df.list$coverage <- cbind(CpG.sum.df.list$coverage, avg=rowMeans(CpG.sum.df.list$coverage[,-1],na.rm=TRUE))
  
  if(length(other.Cs.summary) > 0){
    other.Cs.sum.df.list <- lapply(1:3, function(var){
      tab <- Reduce(function(x,y) full_join(x,y), lapply(other.Cs.summary, function(x){x[[var]]}))
      arrange(tab,id)
    })
    names(other.Cs.sum.df.list) <- c("coverage","freq.meth","errors")
    other.Cs.sum.df.list$coverage <- cbind(other.Cs.sum.df.list$coverage, avg=rowMeans(other.Cs.sum.df.list$coverage[,-1],na.rm=TRUE))
  } else other.Cs.sum.df.list <- NULL
  
  return(list(CpG.sum = CpG.sum.df.list, non.CpG.sum = other.Cs.sum.df.list))
}  
