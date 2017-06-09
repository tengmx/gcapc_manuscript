###### refine peaks using gcapc

library(gcapc) ## version 1.0.8
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
bams <- list.files(".",pattern='.+bam$')
for(bam in bams){
    coverage <- read5endCoverage(bam,chroms=paste0('chr',c(1:22,'X')))
    bdwidth <- bindWidth(cov,range = c(50L, 400L)) 
    gcbias <- gcEffects(cov,bdwidth,plot=T,sampling=c(0.05,5),theta0=1,theta1=10,converge=1e-3)
    name <- strsplit(bam,"\\.")[[1]][1]
    ## macs2
    macs2peaks <- read.table(paste0("macs2/",name,"_peaks.narrowPeak"))
    macs2peaks <- GRanges(macs2peaks$V1,IRanges(macs2peaks$V2+1,macs2peaks$V3),pv=macs2peaks$V8)
    seqlevels(macs2peaks,force=T) <- paste0('chr',c(1:22,'X'))
    seqlengths(macs2peaks) <- seqlengths(Hsapiens)[seqlevels(macs2peaks)]
    refpeaks <- refinePeaks(coverage,gcbias,bdwidth,macs2peaks)
    save(refpeaks,file=paste0("refine/",gsub(".bam","",bam),".macs2.rda"))
    ## hotspot
    hotspotpeaks <- read.table(paste0("hotspot/",name,"/",name,"-final/",name,".hot.bed"))
    hotspotpeaks <- GRanges(hotspotpeaks$V1,IRanges(hotspotpeaks$V2+1,hotspotpeaks$V3),sc=hotspotpeaks$V5)
    seqlevels(hotspotpeaks,force=T) <- paste0('chr',c(1:22,'X'))
    seqlengths(hotspotpeaks) <- seqlengths(Hsapiens)[seqlevels(hotspotpeaks)]
    refpeaks <- refinePeaks(coverage,gcbias,bdwidth,hotspotpeaks)
    save(refpeaks,file=paste0("refine/",gsub(".bam","",bam),".hotspot.rda"))
}
