library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
load('ensembl_m1.rda') ## human CTCF M1 motif by Schmidt et al.
pwm <- round(PWM(M1$int),3)
samples <- c("wgEncodeBroadHistoneHuvecCtcfStdRawDataRep1","wgEncodeBroadHistoneHuvecCtcfStdRawDataRep2",
             "wgEncodeOpenChromChipHuvecCtcfRawDataRep1","wgEncodeOpenChromChipHuvecCtcfRawDataRep2",
             "wgEncodeUwTfbsHuvecCtcfStdRawDataRep1","wgEncodeUwTfbsHuvecCtcfStdRawDataRep2")
library(GenomicRanges)
## loading peak calling results
setwd("~/workspace/chipseq/gcapc/rdas/")
gcapc <- list()
for(i in seq_along(samples)){
    load(paste0("201611292.",samples[i],".peaks.rda"))
    seqlevels(peaks,force=T) <- paste0('chr',c(1:22,'X'))
    peaks <- peaks[peaks$pv<=0.1]
    gcapc[[i]] <- peaks
}
setwd("~/workspace/chipseq/gcapc/spp/0.99/")
spp <- list()
for(i in seq_along(samples)){
    tmp <- read.table(paste0(samples[i],".peak"))
    pr <- GRanges(tmp$V1,IRanges(tmp$V2+1,tmp$V3),sig=tmp$V7,qv=tmp$V9)
    seqlevels(pr,force=T) <- paste0('chr',c(1:22,'X'))
    seqlengths(pr) <- seqlengths(Hsapiens)[seqlevels(pr)]
    pr <- pr[pr$qv>=1]
    spp[[i]] <- trim(pr)
}
# calculate pwm score
pwmsc <- function(peaks,pwm){
    library(BSgenome.Hsapiens.UCSC.hg19)
    library(matrixStats)
    seqs <- getSeq(Hsapiens,peaks)
    possc <- sapply(seqs,function(x) max(mcols(matchPWM(pwm,x,min.score="0%",with.score=TRUE))$score))
    negsc <- sapply(seqs,function(x) max(mcols(matchPWM(reverseComplement(pwm),x,min.score="0%",with.score=TRUE))$score))
    rowMaxs(cbind(possc,negsc))
}
gcapcpwm <- spppwm <- list()
for(i in seq_along(gcapc)){
    #gcapcpwm[[i]] <- pwmsc(gcapc[[i]],pwm)
    spppwm[[i]] <- pwmsc(spp[[i]],pwm)
    cat(i,'\n')
    save(gcapcpwm,spppwm,gcapc,spp,file='~/share/gcapc/20170109.pwmscore.pv0.1.gcapc.spp.rda')
}
# pwm cutoff
layout(matrix(1:6,2,3,byrow=T))
steps <- seq(1000,40000,200)
for(cutoff in seq(0.67,0.92,0.05)){
    for(i in 1:6){
        index1 <- order(gcapc[[i]]$es,decreasing=T)
        index2 <- order(spp[[i]]$sig,decreasing=T)
        prop1 <- sapply(steps,function(x) sum(gcapcpwm[[i]][index1[seq_len(x)]]>=cutoff)/x)
        prop2 <- sapply(steps,function(x) sum(spppwm[[i]][index2[seq_len(x)]]>=cutoff)/x)
        if(i==1) {
            plot(steps,prop1,type='l',main=paste0('pwm cutoff:',cutoff),
                 xlab='ranked peak #',ylab='prop of >cutoff',
                 ylim=c(min(c(prop1,prop2)),max(c(prop1,prop2))))
            lines(steps,prop2,lty=2)
        }else {
            lines(steps,prop1)
            lines(steps,prop2,lty=2)
        }
    }
}
# figure 6B
par(cex.axis=1.2,cex.lab=1.2,font.lab=2,cex.main=1.2,cex=1.1,lwd=2)
colors <- c("#e66101","#999999")
cutoff <- 0.72
steps <- seq(1000,30000,200)
for(i in 1:6){
    index1 <- order(gcapc[[i]]$es,decreasing=T)
    index2 <- order(spp[[i]]$sig,decreasing=T)
    prop1 <- 1-sapply(steps,function(x) sum(gcapcpwm[[i]][index1[seq_len(x)]]>=cutoff)/x)
    prop2 <- 1-sapply(steps,function(x) sum(spppwm[[i]][index2[seq_len(x)]]>=cutoff)/x)
    if(i==1) {
        plot(steps,prop1,type='l',
             xlab='# of ranked peaks by significance ',ylab='False positive rate',
             ylim=c(min(c(prop1,prop2)),max(c(prop1,prop2))),col=colors[1])
        lines(steps,prop2,lty=2,col=colors[2])
    }else {
        lines(steps,prop1,col=colors[1])
        lines(steps,prop2,lty=2,col=colors[2])
    }
}
# numbers
steps <- c(100,500,1000,2000,5000,10000,30000)
for(i in 1:6){
    index1 <- order(gcapc[[i]]$es,decreasing=T)
    index2 <- order(spp[[i]]$sig,decreasing=T)
    prop1 <- sapply(steps,function(x) x-sum(gcapcpwm[[i]][index1[seq_len(x)]]>=cutoff))
    prop2 <- sapply(steps,function(x) x-sum(spppwm[[i]][index2[seq_len(x)]]>=cutoff))
    cat('rep',i,'gcapc',prop1,'spp',prop2,'\n')
}
