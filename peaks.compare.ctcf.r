###### consistency, motif enrichment, refine & sampling proportion  analysis for CTCF

library(matrixStats)
library(GenomicRanges)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(RColorBrewer)

## load called peaks
samples <- c("wwgEncodeUwTfbsHuvecCtcfStdRawDataRep1","wgEncodeUwTfbsHuvecCtcfStdRawDataRep2",
             "wgEncodeOpenChromChipHuvecCtcfRawDataRep1","wgEncodeOpenChromChipHuvecCtcfRawDataRep2",
             "gEncodeBroadHistoneHuvecCtcfStdRawDataRep1","wgEncodeBroadHistoneHuvecCtcfStdRawDataRep2")
labs <- c('Broad','Broad','UTA','UTA','UW','UW')
setwd("gcapc")
# 051: 5%*1; 055: 5%*5; 151: 15%*1; 251: 25%*1; 055T: supervised 5%*5
gcapc051 <- gcapc055 <- gcapc151 <- gcapc251 <- gcapc055T <- list()
for(i in seq_along(samples)){ ## load gcapc peaks
    load(paste0(samples[i],".0.05.1.FALSE.rda"))
    seqlevels(peaks,force=T) <- paste0('chr',c(1:22,'X'))
    seqlengths(peaks) <- seqlengths(Hsapiens)[seqlevels(peaks)]
    peaks <- peaks[order(peaks$es,decreasing=T)][1:min(30000,length(peaks))]
    gcapc051[[i]] <- peaks
    load(paste0(samples[i],".0.05.5.FALSE.rda"))
    seqlevels(peaks,force=T) <- paste0('chr',c(1:22,'X'))
    seqlengths(peaks) <- seqlengths(Hsapiens)[seqlevels(peaks)]
    peaks <- peaks[order(peaks$es,decreasing=T)][1:min(30000,length(peaks))]
    gcapc055[[i]] <- peaks
    load(paste0(samples[i],".0.15.1.FALSE.rda"))
    seqlevels(peaks,force=T) <- paste0('chr',c(1:22,'X'))
    seqlengths(peaks) <- seqlengths(Hsapiens)[seqlevels(peaks)]
    peaks <- peaks[order(peaks$es,decreasing=T)][1:min(30000,length(peaks))]
    gcapc151[[i]] <- peaks
    load(paste0(samples[i],".0.25.1.FALSE.rda"))
    seqlevels(peaks,force=T) <- paste0('chr',c(1:22,'X'))
    seqlengths(peaks) <- seqlengths(Hsapiens)[seqlevels(peaks)]
    peaks <- peaks[order(peaks$es,decreasing=T)][1:min(30000,length(peaks))]
    gcapc251[[i]] <- peaks
    load(paste0(samples[i],".0.05.5.TRUE.rda"))
    seqlevels(peaks,force=T) <- paste0('chr',c(1:22,'X'))
    seqlengths(peaks) <- seqlengths(Hsapiens)[seqlevels(peaks)]
    peaks <- peaks[order(peaks$es,decreasing=T)][1:min(30000,length(peaks))]
    gcapc055T[[i]] <- peaks
}
setwd("spp")
spp <- list()  ## load SPP peaks
for(i in seq_along(samples)){
    peaks <- read.table(paste0(samples[i],".peak"))
    peaks <- GRanges(peaks$V1,IRanges(peaks$V2+1,peaks$V3),sig=peaks$V7)
    seqlevels(peaks,force=T) <- paste0('chr',c(1:22,'X'))
    seqlengths(peaks) <- seqlengths(Hsapiens)[seqlevels(peaks)]
    peaks <- peaks[order(peaks$sig,decreasing=T)][1:min(30000,length(peaks))]
    spp[[i]] <- trim(peaks)
}
setwd("macs2")
macs2 <- list() ## load MACS peaks
for(i in seq_along(samples)){
    peaks <- read.table(paste0(samples[i],"_peaks.narrowPeak"),sep='\t',header=F,stringsAsFactors=F)
    peaks <- GRanges(peaks$V1,IRanges(start=peaks$V2+1,end=peaks$V3),pv=peaks$V8)
    seqlevels(peaks,force=T) <- paste0('chr',c(1:22,'X'))
    seqlengths(peaks) <- seqlengths(Hsapiens)[seqlevels(peaks)]
    peaks <- peaks[order(peaks$pv,decreasing=T)][1:min(30000,length(peaks))]
    macs2[[i]] <- peaks
}
setwd("hotspot")
hotspot <- list() ## load hotspot peaks
for(i in seq_along(samples)){
    peaks <- read.table(paste0(samples[i],"/",samples[i],"-final/",samples[i],".hot.bed"))
    peaks <- GRanges(peaks$V1,IRanges(start=peaks$V2+1,end=peaks$V3),sc=peaks$V5)
    seqlevels(peaks,force=T) <- paste0('chr',c(1:22,'X'))
    seqlengths(peaks) <- seqlengths(Hsapiens)[seqlevels(peaks)]
    peaks <- peaks[order(peaks$sc,decreasing=T)][1:min(30000,length(peaks))]
    hotspot[[i]] <- peaks
}
setwd("refine")
macs2ref <- hotspotref <- list() ## load refined peaks for MACS and hotspot
for(i in seq_along(samples)){
    load(paste0(samples[i],".macs2.rda"))
    seqlevels(refpeaks,force=T) <- paste0('chr',c(1:22,'X'))
    seqlengths(refpeaks) <- seqlengths(Hsapiens)[seqlevels(refpeaks)]
    refpeaks <- refpeaks[order(refpeaks$newes,decreasing=T)][1:min(30000,length(refpeaks))]
    macs2ref[[i]] <- refpeaks
    load(paste0(samples[i],".hotspot.rda"))
    seqlevels(refpeaks,force=T) <- paste0('chr',c(1:22,'X'))
    seqlengths(refpeaks) <- seqlengths(Hsapiens)[seqlevels(refpeaks)]
    refpeaks <- refpeaks[order(refpeaks$newes,decreasing=T)][1:min(30000,length(refpeaks))]
    hotspotref[[i]] <- refpeaks
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
## plot consistency between gcapc/refined and spp/MACS/hotspot
compare2(gcapc055,spp,withinlab=FALSE,main="",labs=labs,legends=c("gcapc","SPP"),ranks=seq(800,30000,50),ylim=c(0,0.9))
compare2(gcapc055,macs2,withinlab=FALSE,main="",labs=labs,legends=c("gcapc","MACS2"),ranks=seq(800,30000,50),ylim=c(0,0.9))
compare2(macs2ref,macs2,ip1=2,ip2=1,withinlab=FALSE,main="",labs=labs,legends=c("MACS2 refined","MACS2"),
         ranks=seq(800,30000,50),ylim=c(0,0.9))
compare2(hotspotref,hotspot,ip1=2,ip2=1,withinlab=FALSE,main="",labs=labs,legends=c("hotspot refined","hotspot"),
         ranks=seq(800,30000,50),ylim=c(0,0.9))

## consistency when using different sampling proportions 
cols <- brewer.pal(7,'Set2')
peakcat(gcapc051[[1]],gcapc151[[1]],col=cols[1],ylim=c(0.95,1),ylab="Proportion in common", xlab="Size of list",
        ranks=seq(800,30000,100))
peakcat(gcapc051[[]],gcapc251[[1]],col=cols[2],add=T,ranks=seq(800,30000,100))
peakcat(gcapc051[[1]],gcapc055[[1]],col=cols[3],add=T,ranks=seq(800,30000,100))
peakcat(gcapc151[[1]],gcapc251[[1]],col=cols[4],add=T,ranks=seq(800,30000,100))
peakcat(gcapc151[[1]],gcapc055[[1]],col=cols[5],add=T,ranks=seq(800,30000,100))
peakcat(gcapc251[[1]],gcapc055[[1]],col=cols[6],add=T,ranks=seq(800,30000,100))
peakcat(gcapc055[[1]],gcapc055T[[1]],col=cols[7],add=T,ranks=seq(800,30000,100))
legend('bottomright',c("5% vs 15%","5% vs 25%","5% vs 5%*5","15% vs 25%","15% vs 5%*5","25% vs 5%*5","5%*5 vs 5%*5 (supervised)"),
       col=cols,lty=1,bty='n')

## motif enrichments calculation
load('JASPAR.rda')
pwm <- round(PWM(CTCF.MA0139.1,prior.params = c(A=0.3, C=0.2, G=0.2, T=0.3)),3)
pwmsc <- function(peaks,pwm){
    seqs <- getSeq(Hsapiens,peaks)
    possc <- sapply(seqs,function(x) max(mcols(matchPWM(pwm,x,min.score="0%",with.score=TRUE))$score))
    negsc <- sapply(seqs,function(x) max(mcols(matchPWM(reverseComplement(pwm),x,min.score="0%",with.score=TRUE))$score))
    rowMaxs(cbind(possc,negsc))
}
gcapc <- gcapc055
gcapcpwm <- spppwm <- macs2pwm <- list()
for(i in seq_along(gcapc)){
    gcapcpwm[[i]] <- pwmsc(gcapc[[i]],pwm)
    spppwm[[i]] <- pwmsc(spp[[i]],pwm)
    macs2pwm[[i]] <- pwmsc(macs2[[i]],pwm)
    cat(i,'\n')
}

## motif enrichment performance independent of pwm cutoff: gcapc vs SPP
layout(matrix(1:6,2,3,byrow=T))
steps <- seq(1000,30000,200)
colors <- c("#e66101","#999999")
for(cutoff in seq(0.67,0.92,0.05)){
    for(i in 1:6){
        index1 <- order(gcapc[[i]]$es,decreasing=T)
        index2 <- order(spp[[i]]$sig,decreasing=T)
        prop1 <- 1-sapply(steps,function(x) sum(gcapcpwm[[i]][index1[seq_len(x)]]>=cutoff)/x)
        prop2 <- 1-sapply(steps,function(x) sum(spppwm[[i]][index2[seq_len(x)]]>=cutoff)/x)
        if(i==1) {
            plot(steps,prop1,type='l',main=paste0('pwm cutoff:',cutoff),
                 xlab='ranked peak #',ylab='False positive rate',col=colors[1],
                 ylim=c(0,max(c(prop1,prop2))))
            lines(steps,prop2,lty=2,col=colors[2])
        }else {
            lines(steps,prop1,col=colors[1])
            lines(steps,prop2,lty=2,col=colors[2])
        }
    }
}
legend('bottomright',c('gcapc','SPP'),col=colors,lty=1:2,bty='n')

## motif enrichmen comparison between gcapc and MACS
colors <- c("#e66101","#999999")
cutoff <- 0.72
steps <- seq(1000,30000,200)
for(i in 1:6){
    index1 <- order(gcapc[[i]]$es,decreasing=T)
    index2 <- order(macs2[[i]]$pv,decreasing=T)
    prop1 <- 1-sapply(steps,function(x) sum(gcapcpwm[[i]][index1[seq_len(x)]]>=cutoff)/x)
    prop2 <- 1-sapply(steps,function(x) sum(macs2pwm[[i]][index2[seq_len(x)]]>=cutoff)/x)
    if(i==1) {
        plot(steps,prop1,type='l',
             xlab='# of ranked peaks by significance ',ylab='False positive rate',
             ylim=c(0,0.018),col=colors[1])
        lines(steps,prop2,lty=2,col=colors[2])
    }else {
        lines(steps,prop1,col=colors[1])
        lines(steps,prop2,lty=2,col=colors[2])
    }
}
legend('bottomright',c('gcapc','MACS2'),col=colors,lty=1:2,bty='n')

##  total false positives
steps <- c(100,500,1000,2000,5000,10000,30000)
cutoff <- 0.72
for(i in 1:6){
    index1 <- order(gcapc[[i]]$es,decreasing=T)
    index2 <- order(spp[[i]]$sig,decreasing=T)
    prop1 <- sapply(steps,function(x) x-sum(gcapcpwm[[i]][index1[seq_len(x)]]>=cutoff))
    prop2 <- sapply(steps,function(x) x-sum(spppwm[[i]][index2[seq_len(x)]]>=cutoff))
    cat('rep',i,'gcapc',prop1,'spp',prop2,'\n')
}
