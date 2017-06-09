###### consistency analysis for EP300

library(matrixStats)
library(GenomicRanges)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)

## load called peaks
bams <- c('wgEncodeHaibTfbsHepg2P300V0416101RawDataRep1.bam','wgEncodeHaibTfbsHepg2P300V0416101RawDataRep2.bam',
          'wgEncodeSydhTfbsHepg2P300sc582IggrabRawDataRep1.bam','wgEncodeSydhTfbsHepg2P300sc582IggrabRawDataRep2.bam')
samples <- gsub('.bam','',bams)
labs <- c('HAIB','HAIB','SYDH','SYDH')
num <- 10000
setwd("gcapc")
gcapc <- list()
for(i in seq_along(samples)){
    load(paste0(samples[i],".rda"))
    seqlevels(peaks,force=T) <- paste0('chr',c(1:22,'X'))
    seqlengths(peaks) <- seqlengths(Hsapiens)[seqlevels(peaks)]
    peaks <- peaks[order(peaks$es,decreasing=T)][1:min(num,length(peaks))]
    gcapc[[i]] <- peaks
}
setwd("spp")
spp <- list()
for(i in seq_along(samples)){
    peaks <- read.table(paste0(samples[i],".peak"))
    peaks <- GRanges(peaks$V1,IRanges(peaks$V2+1,peaks$V3),sig=peaks$V7)
    seqlevels(peaks,force=T) <- paste0('chr',c(1:22,'X'))
    seqlengths(peaks) <- seqlengths(Hsapiens)[seqlevels(peaks)]
    peaks <- peaks[order(peaks$sig,decreasing=T)][1:min(num,length(peaks))]
    spp[[i]] <- peaks
}
peakcat <- function(x,y,ranks=seq(200,20000,50),ix=1,iy=1,add=FALSE,...){
    fo <- findOverlaps(x,y)
    xfo <- queryHits(fo)
    yfo <- subjectHits(fo)
    xo <- order(mcols(x)[,ix],decreasing=TRUE)
    yo <- order(mcols(y)[,iy],decreasing=TRUE)
    xycat <- sapply(ranks,function(t){
        idx <- seq_len(t)
        # minor bug, trival
        sum(xfo %in% xo[idx] & yfo %in% yo[idx])/t
    })
    if(add) lines(ranks,xycat,...)
    else plot(ranks,xycat,type='l',...)
}
compare2 <- function(p1,p2,ip1=1,ip2=1,legends=c("gcapc","XXXX"),
                     ranks=seq(800,20000,50),labs,withinlab=TRUE,...){
    library(RColorBrewer)
    colors <- c("#e66101","#999999")
    if(length(p1)!=length(p2)) stop("length error!")
    p1med <- sapply(p1,function(x) median(width(x)))
    p2med <- sapply(p2,function(x) median(width(x)))
    for(i in seq_along(p1)){
        tmp <- p2med[i]/p1med[i]
        resized <- round(width(p2[[i]])/tmp)
        shifted <- round((width(p2[[i]])-resized)/2)
        p2[[i]] <- shift(resize(p2[[i]],resized),shifted)
    }
    plot(NULL, xlim=c(min(ranks),max(ranks)),
         ylab="Proportion in common", xlab="Size of list",...)
    for( i in seq_len(length(p1)-1)){
        for(j in (i+1):length(p1)){
            cond <- labs[i]==labs[j]
            if(!withinlab){
                cond <- !cond
            }
            if(cond){
                peakcat(x=p1[[i]],y=p1[[j]],ranks=ranks[ranks <= min(length(p1[[i]]),length(p1[[j]]))],
                        ix=ip1,iy=ip1,add=T,col=colors[1],lty=1)
                peakcat(x=p2[[i]],y=p2[[j]],ranks=ranks[ranks <= min(length(p2[[i]]),length(p2[[j]]))],
                        ix=ip2,iy=ip2,add=T,col=colors[2],lty=2)
                cat(i,j,'\n')
            }
        }
    }
    legend('bottomright',legends,col=colors,lty=1:2,bty='n')
}
compare2(gcapc,spp,withinlab=FALSE,main="EP300 HepG2",labs=labs,legends=c("gcapc","SPP"),ranks=seq(800,10000,50),ylim=c(0,0.6))

## gc content of called peaks
gccont.human <- function(grg){
    grgseq <- getSeq(Hsapiens,grg)
    wsf <- alphabetFrequency(grgseq, baseOnly=T, as.prob=T)
    rowSums(wsf[,c("C","G")])/rowSums(wsf)
}
gcapcgc <- sppgc <- list()
for(i in 1:4){
    gcapcgc[[i]] <- gccont.human(gcapc[[i]])
    sppgc[[i]] <- gccont.human(spp[[i]])
}
boxplot(gcapcgc,main="",names=labs,las=2,ylab='GC content of gcapc peaks',pch=20,cex=0.3,ylim=c(0.2,1))
boxplot(sppgc,main="",names=labs,las=2,ylab='GC content of SPP peaks',pch=20,cex=0.3,ylim=c(0.2,1))
