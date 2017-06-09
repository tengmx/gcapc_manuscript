###### consistency and motif analysis for GATA2

library(matrixStats)
library(GenomicRanges)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)

## load called peaks
samples <- c('wgEncodeHaibTfbsK562Gata2sc267Pcr1xRawDataRep1','wgEncodeHaibTfbsK562Gata2sc267Pcr1xRawDataRep2',
             'wgEncodeSydhTfbsK562Gata2RawDataRep1','wgEncodeSydhTfbsK562Gata2RawDataRep2')
labs <- c('HAIB','HAIB','SYDH','SYDH')
num <- 10000
setwd("gcapc")
gcapc <- list()
for(i in seq_along(samples)){
    load(paste0(samples[i],".rda"))
    seqlevels(peaks,force=T) <- paste0('chr',c(1:22,'X'))
    seqlengths(peaks) <- seqlengths(Hsapiens)[seqlevels(peaks)]
    peaks <- peaks[order(peaks$es,decreasing=T)][1:min(10000,length(peaks))]
    gcapc[[i]] <- peaks
}
setwd("spp")
spp <- list()
for(i in seq_along(samples)){
    tmp <- read.table(paste0(samples[i],".peak"))
    pr <- GRanges(tmp$V1,IRanges(tmp$V2+1,tmp$V3),sig=tmp$V7)
    seqlevels(pr,force=T) <- paste0('chr',c(1:22,'X'))
    seqlengths(pr) <- seqlengths(Hsapiens)[seqlevels(pr)]
    pr <- pr[order(pr$sig,decreasing=T)][1:min(10000,length(peaks))]
    spp[[i]] <- pr
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
compare2(gcapc,spp,withinlab=FALSE,main="GATA2 K562",labs=labs,legends=c("gcapc","SPP"),ranks=seq(800,10000,50),ylim=c(0,0.7))

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

## motif enrichment
load('JASPAR.rda')
pwmsc <- function(peaks,pwm){
    seqs <- getSeq(Hsapiens,peaks)
    possc <- sapply(seqs,function(x) max(mcols(matchPWM(pwm,x,min.score="0%",with.score=TRUE))$score))
    negsc <- sapply(seqs,function(x) max(mcols(matchPWM(reverseComplement(pwm),x,min.score="0%",with.score=TRUE))$score))
    rowMaxs(cbind(possc,negsc))
}
pwm <- round(PWM(GATA2.MA0036.2,prior.params = c(A=0.3, C=0.2, G=0.2, T=0.3)),3)
gcapcpwm <- spppwm <- list()
for(i in seq_along(gcapc)){
    gcapcpwm[[i]] <- pwmsc(gcapc[[i]],pwm)
    spppwm[[i]] <- pwmsc(spp[[i]],pwm)
    cat(i,'\n')
}
colors <- c("#e66101","#999999")
cutoff <- 0.72
steps <- seq(1000,10000,200)
for(i in 1:4){
    index1 <- order(gcapc[[i]]$es,decreasing=T)
    index2 <- order(spp[[i]]$sig,decreasing=T)
    prop1 <- 1-sapply(steps,function(x) sum(gcapcpwm[[i]][index1[seq_len(x)]]>=cutoff)/x)
    prop2 <- 1-sapply(steps,function(x) sum(spppwm[[i]][index2[seq_len(x)]]>=cutoff)/x)
    if(i==1) {
        plot(steps,prop1,type='l',
             xlab='# of ranked peaks by significance ',ylab='False positive rate',
             ylim=c(0,0.12),col=colors[1],lty=2)
        lines(steps,prop2,lty=2,col=colors[2])
    }else {
        lines(steps,prop1,col=colors[1],lty=2)
        lines(steps,prop2,lty=2,col=colors[2])
    }
}
legend('topleft',c('gcapc','SPP'),col=colors,lty=2,bty='n')
