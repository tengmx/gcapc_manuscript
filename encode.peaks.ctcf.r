###### CTCF peaks reported by ENCODE

library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(VennDiagram)

## load peaks downloaded from ENCODE portal
peakfiles <- list.files(".",pattern="^Huvec.+bed$")
peaks <- GRangesList()
for(i in seq_along(peakfiles)){
    tmp <- read.table(peakfiles[i],sep='\t',header=F,stringsAsFactors=F)
    tmp <- makeGRangesFromDataFrame(tmp,seqnames.field="V1",start.field="V2",end.field="V3",
                                    strand.field="V6",starts.in.df.are.0based=TRUE,keep.extra.columns=TRUE)
    peaks[[i]] <- tmp
    cat(i, '\n')
}

## discordance
cellpeaks <- reduce(do.call(c,peaks))
fo <- lapply(peaks,function(x) findOverlaps(x,cellpeaks))
foidxlst <- lapply(fo,function(x) unique(subjectHits(x)))
foidx <- table(unlist(foidxlst))

## GC content
gccont.human <- function(grg){
    grgseq <- getSeq(Hsapiens,grg)
    wsf <- alphabetFrequency(grgseq, baseOnly=T, as.prob=T)
    rowSums(wsf[,c("C","G")])/rowSums(wsf)
}
gc <- gccont.human(cellpeaks)
gcs <- lapply(1:3,function(i) gc[intersect(which(foidx==1),foidxlst[[i]])])
boxplot(c(gcs,list(gc[foidx==3])),ylab='GC content',pch=20,cex=0.3,ylim=c(0.2,0.86))

## venn diagram: peaks are mergerd to accounts for multi-to-single overlapping
gcs1 <- sapply(1:3,function(i) length(intersect(which(foidx==1),foidxlst[[i]])))
gcs2 <- sapply(1:3,function(i) length(intersect(intersect(which(foidx==2),foidxlst[[i]]),foidxlst[[i%%3+1]])))
gcs3 <- length(which(foidx==3))
grid.newpage()
draw.triple.venn(area1 = gcs3+gcs1[3]+gcs2[2]+gcs2[3],
                 area2 = gcs3+gcs1[2]+gcs2[1]+gcs2[2],
                 area3 = gcs3+gcs1[1]+gcs2[1]+gcs2[3], n12 = gcs3+gcs2[2], n23 = gcs3+gcs2[1], n13 = gcs3+gcs2[3],
                 n123 = gcs3, category = c("UW", "UTA", "Broad"), lwd=2,
                 fill = NULL)

## proportion of only by one lab
sum(gcs1)/sum(gcs1,gcs2,gcs3)
