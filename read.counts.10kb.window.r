###### counting reads for 10kb genomic windows

library(Rsubread)
library(GenomicRanges)
library(Rsamtools)
library(BSgenome.Hsapiens.UCSC.hg19)

## samples for trancription factor
samples <- list.files(".",pattern="^wgEncode.+.bam$",full=T)
binsize <- 10000

## 10kb bins
seqs <- sortSeqlevels(seqinfo(BamFileList(samples)))
seqs <- seqlengths(seqs)[!grepl("_", seqnames(seqs))]
starts <- lapply(seqs, function(i) seq(1, i, binsize))
ends <- lapply(seqs, function(i) {
    if(i%%binsize == 0) seq(binsize, i, binsize)
    else c(seq(binsize, i, binsize), i)
})
chrs <- rep(names(seqs), times=sapply(starts, length))
regions <- GRanges(chrs, IRanges(start=unlist(starts),
                                 end=unlist(ends)))
seqlengths(regions) <- seqs[seqlevels(regions)]

## count reads by overlapping
anno <- data.frame(GeneID=seq_len(length(regions)),
                   Chr=seqnames(regions), Start=start(regions),
                   End=end(regions), Strand=strand(regions))
count <- featureCounts(files=samples, annot.ext=anno,
                       read2pos='5', ignoreDup=TRUE)
## GC content calculation
seqs <- Views(Hsapiens,regions)
gcpos <- alphabetFrequency(seqs)
gc <- rowSums(gcpos[,c('C','G')])/binsize
gc[gc<0.28 | gc>0.67] <- NA

## plot
for(i in seq_len(ncol(count$counts))){
    plot(gc,count$counts[,i]+1,log='y',col=rgb(0,0,0,alpha=0.05),pch=20,ylim=c(1,2000),
         xlab='GC content',ylab="Read counts + 1")
    cat(i,'\n')
}
