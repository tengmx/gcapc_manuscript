samples <- c("wgEncodeBroadHistoneHuvecCtcfStdRawDataRep1","wgEncodeBroadHistoneHuvecCtcfStdRawDataRep2",
             "wgEncodeOpenChromChipHuvecCtcfRawDataRep1","wgEncodeOpenChromChipHuvecCtcfRawDataRep2",
             "wgEncodeUwTfbsHuvecCtcfStdRawDataRep1","wgEncodeUwTfbsHuvecCtcfStdRawDataRep2")
library(GenomicRanges)
setwd("~/workspace/chipseq/gcapc/rdas/")
gcapc <- list()
for(i in seq_along(samples)){
    load(paste0("201611292.",samples[i],".peaks.rda"))
    seqlevels(peaks,force=T) <- paste0('chr',c(1:22,'X'))
    gcapc[[i]] <- peaks
}
setwd("~/workspace/chipseq/gcapc/spp/0.99/")
spp <- list()
for(i in seq_along(samples)){
    tmp <- read.table(paste0(samples[i],".peak"))
    pr <- GRanges(tmp$V1,IRanges(tmp$V2+1,tmp$V3),sig=tmp$V7)
    seqlevels(pr,force=T) <- paste0('chr',c(1:22,'X'))
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
                     ranks=seq(800,20000,50),withinlab=TRUE,...){
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
            cond <- !(i %in% seq(1,length(p1),2)) || i+1 != j
            if(withinlab){
                cond <- !cond
            }
            if(cond){
                peakcat(x=p1[[i]],y=p1[[j]],ranks=ranks,ix=ip1,iy=ip1,add=T,col=colors[1],lty=1)
                peakcat(x=p2[[i]],y=p2[[j]],ranks=ranks,ix=ip2,iy=ip2,add=T,col=colors[2],lty=2)
                cat(i,j,'\n')
            }
        }
    }
    legend('bottomright',legends,col=colors,lty=1:2,bty='n')
}
par(cex.axis=1.2,cex.lab=1.2,font.lab=2,cex.main=1.2,cex=1.1,lwd=2)
compare2(gcapc,spp,withinlab=FALSE,main="",legends=c("gcapc","SPP"),ranks=seq(800,30000,50),ylim=c(0,0.9))



