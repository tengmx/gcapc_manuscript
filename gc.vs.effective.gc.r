###### comparing GC content and Effective GC content for bins

library(BSgenome.Hsapiens.UCSC.hg19)

## GC content with different bin size
gcct <- GRangesList()
binsizes <- c(200L,300L,400L,500L)
for(i in seq_along(binsizes)){
    binsize <- binsizes[i]
    seqs <- seqlengths(Hsapiens)[1:24]
    starts <- lapply(seqs, function(i) seq(1, i, binsize))
    ends <- lapply(seqs, function(i) {
        if(i%%binsize == 0) seq(binsize, i, binsize)
        else c(seq(binsize, i, binsize), i)
    })
    chrs <- rep(names(seqs), times=sapply(starts, length))
    regions <- GRanges(chrs, IRanges(start=unlist(starts),
                                     end=unlist(ends)))
    seqlengths(regions) <- seqs[seqlevels(regions)]
    seqs <- getSeq(Hsapiens,regions)
    gcpos <- alphabetFrequency(seqs)
    gc <- rowSums(gcpos[,c('C','G')])/binsize
    mcols(regions)$gc <- gc
    gcct[[i]] <- regions
    cat(i,'\n')
}
for(i in seq_along(gcct)){
    mcols(gcct[[i]])$gc[mcols(gcct[[i]])$gc==0] <- NA
}
rgidx <- seq_along(gcct[[length(gcct)]])
gcdat <- matrix(0,length(rgidx),length(gcct))
for(i in seq_len(length(gcct)-1)){
    fo <- findOverlaps(gcct[[i]],gcct[[length(gcct)]])
    gcdat[,i] <- mcols(gcct[[i]])$gc[queryHits(fo)[match(rgidx,subjectHits(fo))]]
}
gcdat[,length(gcct)] <- mcols(gcct[[length(gcct)]])$gc

## Effective GC content for 200bp bins with different flanking regions  
bdwidth <- 200
seqs <- seqlengths(Hsapiens)[1:24]
seqs <- floor(seqs/bdwidth-2)*bdwidth
starts <- lapply(seqs, function(i) seq(1+bdwidth*2, i, bdwidth))
ends <- lapply(seqs, function(i) seq(bdwidth*3, i, bdwidth))
chrs <- rep(names(seqs), times=sapply(starts, length))
samp <- 0.05
sampidx <- sort(sample.int(length(chrs),ceiling(length(chrs)*samp)))
region <- GRanges(chrs[sampidx], IRanges(start=unlist(starts)[sampidx],
                                         end=unlist(ends)[sampidx]))
gccont0 <- function(region,flank=200L){
    cat("......... flanking:",flank,'\n')
    rwidth <- width(region[1])
    nr <- shift(resize(region,rwidth + flank*2),-flank)
    seqs <- getSeq(Hsapiens,nr)
    gcpos <- gregexpr('(G)|(C)',seqs)
    weight <- c(seq_len(flank),rep(flank+1,rwidth),rev(seq_len(flank)))
    gcnumw <- sapply(gcpos,function(x) sum(weight[x]))
    gcnumw[gcnumw==(sum(weight)-1)] <- 0
    round(gcnumw / (flank+rwidth) / (flank+1),3)
}
flanks <- c(0,50,100,150)
gcf <- matrix(0,length(region),length(flanks))
for(i in seq_along(flanks)){
    gcf[,i] <- gccont0(region,flank=flanks[i])
    cat(i,'\n')
}
gcf[gcf==0] <- NA

## plot
s1 <- sample(seq_len(nrow(gcdat)),100000)
s2 <- sample(seq_len(nrow(gcf)),100000)
plot(gcdat[s1,2],gcdat[s1,3],col=rgb(0,0,0,alpha=0.2),pch=20,ylim=c(0.2,0.8),xlim=c(0.2,0.8),
     xlab='300bp fragments',ylab="400bp fragments",main="GC Content")
plot(gcf[s2,2],gcf[s2,3],col=rgb(0,0,0,alpha=0.2),pch=20,ylim=c(0.2,0.8),xlim=c(0.2,0.8),
     xlab='300bp region (200bp in + 50*2 out)',ylab="400bp region (200bp in + 100*2 out)",
     main="Effective GC Content")
