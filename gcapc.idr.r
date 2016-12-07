setwd('~/Desktop/gcapc/20161129_4_2/')
library(GenomicRanges)
pfs <- list.files('.',pattern='^201611292.+rda$')
plst <- GRangesList()
for(i in seq_along(pfs)){
  load(pfs[i])
  plst[[i]] <- peaks
  tmp <- as.data.frame(plst[[i]])
  tmp <- data.frame(tmp[,1:3],".",".",tmp[,5:7])
  write.table(tmp,paste0(pfs[i],'.peak'),row.names=F,col.names=F,sep='\t',quote=F)
}
head(as.data.frame(plst[[i]]))

######### using idr in commond line as follows, then generate plots
'
lab=201611292.wgEncodeBroadHistoneHuvecCtcfStdRawData
lab=201611292.wgEncodeOpenChromChipHuvecCtcfRawData
lab=201611292.wgEncodeUwTfbsHuvecCtcfStdRawData
idr --samples $lab\Rep1.peaks.rda.peak \
              $lab\Rep2.peaks.rda.peak \
--input-file-type bed \
--output-file $lab.idr \
--idr-threshold 0.02 \
--plot \
--rank 7
'

pkf <- c('201611292.wgEncodeBroadHistoneHuvecCtcfStdRawData.idr',
         '201611292.wgEncodeOpenChromChipHuvecCtcfRawData.idr',
         '201611292.wgEncodeUwTfbsHuvecCtcfStdRawData.idr')
library(GenomicRanges)
cellpeaklist <- GRangesList()
for(i in seq_along(pkf)){
  tmp <- read.table(pkf[i],sep='\t',header=F,stringsAsFactors=F)
  cellpeaklist[[i]] <- makeGRangesFromDataFrame(tmp,seqnames.field="V1",start.field="V2",end.field="V3",
                                  strand.field="V6",starts.in.df.are.0based=TRUE,keep.extra.columns=FALSE)
  cat(i, '\n')
}
cellpeaks <- reduce(do.call(c,cellpeaklist),min.gapwidth=0)
fo <- lapply(cellpeaklist,function(x) findOverlaps(x,cellpeaks,maxgap = 0,minoverlap = 1))
foidxlst <- lapply(fo,function(x) unique(subjectHits(x)))
foidx <- table(unlist(foidxlst))
densplot <- function(x,adjust=1.2,colors,lty,legend,...){
  for(i in seq_along(x)){
    if(i==1) plot(density(x[[i]],adjust=adjust),col=colors[i],lty=lty[i],...)
    else lines(density(x[[i]],adjust=adjust),col=colors[i],lty=lty[i])
  }
  legend('topright',legend=legend,lty=lty,col=colors,bty='n')
}
gccont.human <- function(grg){
  library(BSgenome.Hsapiens.UCSC.hg19)
  grgseq <- getSeq(Hsapiens,grg)
  wsf <- alphabetFrequency(grgseq, baseOnly=T, as.prob=T)
  rowSums(wsf[,c("C","G")])/rowSums(wsf)
}
gc <- gccont.human(cellpeaks)
library(RColorBrewer)
lab = c('broad','uta','uw')
gcs <- lapply(1:3,function(i) gc[intersect(which(foidx==1),foidxlst[[i]])])
densplot(c(gcs,list(gc[foidx==3])),colors=c(brewer.pal(3,"Set2"),"black"),lty=c(rep(1,3),2),
         xlab='GC content',ylim=c(0,5.2),
         legend=c(paste("Only in",lab),paste("In all labs")))
pdf("~/Desktop/20161130.figure2b.pdf",4.5,6)
par(cex.lab=1.2,cex.axis=1.2,cex.main=1.2,cex=1.2,lwd=2)
boxplot(c(gcs,list(gc[foidx==3]))[c(3,1,2,4)],names=c(paste("Only in",lab),paste("In all labs"))[c(3,1,2,4)],
        ylab='GC content',pch=20,cex=0.3,col=c(brewer.pal(4,"Dark2")[c(4,1,3)],'white'),ylim=c(0.2,0.86))
dev.off()
a1 = sapply(gcs,length)
gcs <- lapply(c(2,3,1),function(i) gc[intersect(which(foidx==2),intersect(foidxlst[[i]],foidxlst[[i%%3+1]]))])
densplot(c(gcs,list(gc[foidx==3])),colors=c(brewer.pal(3,"Set2"),"black"),lty=c(rep(1,3),2),
         xlab='GC content',ylim=c(0,4.8),
         legend=c(paste("Not in",lab),paste("In all labs")))
boxplot(c(gcs,list(gc[foidx==3])),names=c(paste("Not in",lab),paste("In all labs")))
a2 = sapply(gcs,length)
a3 = length(gc[foidx==3])
a3 / (sum(a1)+sum(a2)+a3)
sum(a1)/(sum(a1)+sum(a2)+a3)

library(VennDiagram)
pdf("~/Desktop/20161130.supp.figure2b.pdf",4,4)
par(font.axis=2,font.lab=2,cex.lab=1.2,cex.main=1.2,cex.axis=1.2)
grid.newpage()
draw.triple.venn(area1 = a1[3]+a2[1]+a2[2]+a3, area2 = a1[2]+a2[1]+a2[3]+a3, area3 = a1[1]+a2[2]+a2[3]+a3, 
                 n12 = a2[1]+a3, n23 = a2[3]+a3, n13 = a2[2]+a3,
                 n123 = a3, category = c("UW", "UTA", "Broad"), lwd=2,
                 col = brewer.pal(4,"Dark2")[c(4,3,1)], fill = NULL,print.mode='raw')
dev.off()
