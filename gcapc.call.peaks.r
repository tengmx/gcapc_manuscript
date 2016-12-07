library(gcapc)
reads <- "wgEncodeBroadHistoneHuvecCtcfStdRawDataRep1.bam"
cov <- read5endCoverage(reads,chroms=paste0('chr',c(1:22,'X','M')))
bdwidth <- bindWidth(cov,range = c(50L, 400L))
gcbias <- gcEffects(cov,bdwidth,plot=FALSE,converge=5e-5)
peaks <- gcapcPeaks(cov,gcbias,bdwidth,plot=FALSE,pv=0.99)
save(peaks,file="~/workspace/chipseq/gcapc/rdas/201611292.wgEncodeBroadHistoneHuvecCtcfStdRawDataRep1.peaks.rda")
