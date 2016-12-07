marker <- 'Ctcf'
peakfiles <- list.files(paste0("links/peaks/",marker),pattern="bed$")
meta <- t(as.data.frame(strsplit(peakfiles,"\\."),row.names=c("cell","marker","lab","encodeid","type")))
row.names(meta) <- NULL
library(GenomicRanges)
peaks <- GRangesList()
for(i in seq_along(peakfiles)){
    tmp <- read.table(paste0("links/peaks/",marker,"/",peakfiles[i]),sep='\t',header=F,stringsAsFactors=F)
    tmp <- makeGRangesFromDataFrame(tmp,seqnames.field="V1",start.field="V2",end.field="V3",
                                           strand.field="V6",starts.in.df.are.0based=TRUE,keep.extra.columns=TRUE)
    peaks[[i]] <- tmp
    cat(i, '\n')
}

###### discordance
library(RColorBrewer)
cell <- "Huvec"
lab <- meta[meta[,1]==cell,3]
cellpeaklist <- peaks[meta[,1]==cell]
cellpeaks <- reduce(do.call(c,cellpeaklist))
fo <- lapply(cellpeaklist,function(x) findOverlaps(x,cellpeaks))
foidxlst <- lapply(fo,function(x) unique(subjectHits(x)))
foidx <- table(unlist(foidxlst))
densplot <- function(x,adjust=1.2,colors,lty,legend,...){
    for(i in seq_along(x)){
        if(i==1) plot(density(x[[i]],adjust=adjust),col=colors[i],lty=lty[i],...)
        else lines(density(x[[i]],adjust=adjust),col=colors[i],lty=lty[i])
    }
    legend('topright',legend=legend,lty=lty,col=colors,bty='n')
}
### GC content
gccont.human <- function(grg){
    library(BSgenome.Hsapiens.UCSC.hg19)
    grgseq <- getSeq(Hsapiens,grg)
    wsf <- alphabetFrequency(grgseq, baseOnly=T, as.prob=T)
    rowSums(wsf[,c("C","G")])/rowSums(wsf)
}
gc <- gccont.human(cellpeaks)
par(cex.lab=1.2,cex.axis=1.2,cex.main=1.2,cex=1.2,lwd=2)
gcs <- lapply(1:3,function(i) gc[intersect(which(foidx==1),foidxlst[[i]])])
boxplot(c(gcs,list(gc[foidx==3]))[c(3,1,2,4)],names=c(paste("Only in",lab),paste("In all labs"))[c(3,1,2,4)],
        ylab='GC content',pch=20,cex=0.3,col=c(brewer.pal(4,"Dark2")[c(4,1,3)],'white'),ylim=c(0.2,0.86))

### venn diagram: peaks are mergerd to accounts for multi-to-single overlapping
### thus the sums are not the exact original peak numbers
library(VennDiagram)
par(font.axis=2,font.lab=2,cex.lab=1.2,cex.main=1.2,cex.axis=1.2)
grid.newpage()
draw.triple.venn(area1 = 37847, area2 = 43945, area3 = 36995, n12 = 1760+31212, n23 = 2384+31212, n13 = 2386+31212,
                 n123 = 31212, category = c("UW", "UTA", "Broad"), lwd=2,
                 col = brewer.pal(4,"Dark2")[c(4,3,1)], fill = NULL)

