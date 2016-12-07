###### plot figures

### example: trace
library(GenomicRanges)
load("rdas/Ctcf/v3_sitereadcount/20160314.sites.raw38.q10.ext214.multic.nodup.rda")
load(paste0("rdas/20160314.hg19.gc.150bpsites.400fragsize.rda"))
load("rdas/Ctcf/v3_sitereadcount/20160329.sites.huvec.scale.gccorrect.smallset.rda")
gainidx <- c(19470,13852)
source("~/base/genomic.r")
setwd("~/tmp/Huvec/raw/")
bwf <- list.files(pattern="bw$")
window <- 100
cellmarkidx <- grepl(paste0("CTCF\\(","HUVEC","\\)"), mcols(sites)$detail)
sig <- rawsig[cellmarkidx,meta[,3]=="Huvec"]
gc <- sitegc[cellmarkidx]
met <- meta[meta[,3]=="Huvec",]
sit <- sites[cellmarkidx]
sizefac <- apply(sig,2,quantile,0.75)
par(cex.axis=1.2,cex.lab=1.2,cex.main=1.2,cex=1.1,lwd=2)
layout(matrix(1:2,1,2))
gccontent <- c('Low','High')
for(i in seq_along(gainidx)){
    region <- sit[gainidx[i]]
    start(region) <- start(region) - 600
    end(region) <- end(region) + 600 + window
    cov <- region.coverage.from.bigwigs(bwf,region)
    scalecov <- RleList()
    for(j in seq_len(length(cov))){
        cov[[j]][is.na(cov[[j]])] <- 0
        scalecov[[j]] <- cov[[j]] / sizefac[j] * median(sizefac)
    }
    covplot(scalecov,window,region,main=paste0(gccontent[i]," GC Content"),log=F)
    cat(gainidx[i],'\n')
}
covplot <- function(cov,window,region,main="",log=TRUE){
    library(RColorBrewer)
    colors <- brewer.pal(4,"Dark2")[c(1,3,4)]
    if(log)
        rlelist.track.view.line.log(rlelist.window.smooth(cov,window),
                                    ylim=c(0,40),xaxt='n',
                                    col=rep(colors,each=2),axes=T,
                                    xlab=seqnames(region),ylab="coverage",main=main)
    else
        rlelist.track.view.line(rlelist.window.smooth(cov,window),
                                ylim=c(0,40),xaxt='n',
                                col=rep(colors,each=2),
                                xlab=seqnames(region),ylab="coverage",main=main)
    xaxisidx <- which(start(region):end(region) %% 500 ==0)
    axis(1,at=xaxisidx,labels=xaxisidx + start(region) - 1)
    legend("topleft",c("Broad","UTA","UW"),col=colors,lty=1,bty='n')
}


### region: signal ~ gc  ## need information from previous plot
library(Rsamtools)
library(rtracklayer)
library(Rsubread)
library(BSgenome.Hsapiens.UCSC.hg19)
bams <- c('wgEncodeBroadHistoneHuvecCtcfStdRawDataRep1.bam','wgEncodeBroadHistoneHuvecCtcfStdRawDataRep2.bam',
          'wgEncodeOpenChromChipHuvecCtcfRawDataRep1.bam','wgEncodeOpenChromChipHuvecCtcfRawDataRep2.bam',
          'wgEncodeUwTfbsHuvecCtcfStdRawDataRep1.bam','wgEncodeUwTfbsHuvecCtcfStdRawDataRep2.bam')
source("~/workspace/chipseq/compdiff/ComplexDiff/R/regionCounts.R")
gainidx <- c(19470,13852)
examples <- sit[gainidx]
rc <- regionCounts(bams,binsize=10000L)
seqs <- Views(Hsapiens,rc$regions)
gcpos <- alphabetFrequency(seqs)
gc <- rowSums(gcpos[,c('C','G')])/10000L
gc[gc<0.28 | gc>0.67] <- NA
fo <- findOverlaps(examples,rc$regions)
layout(matrix(1:6,2,3))
maintext <- rep(c('Broad','UTA','UW'),each=2)
par(cex.axis=1.8,cex.lab=1.8,cex.main=1.8,cex=1.1,lwd=2)
library(RColorBrewer)
for(i in c(5,6,3,4,1,2)){
    plot(gc,rc$count[,i]+1,log='y',col=rgb(0,0,0,alpha=0.05),pch=20,ylim=c(1,1000),
         xlab='',ylab="",main="",xaxt='n',yaxt='n')
    points(gc[subjectHits(fo)],rc$count[subjectHits(fo),i]+1,pch=c('G','H'),col=brewer.pal(5,"RdYlBu")[1],cex=2,font=2)
}
dev.off()
par(cex.axis=1.2,cex.lab=1.2,cex.main=1.2,cex=1.1,lwd=2)
plot(NULL,ylim=c(1,1000),xlim=c(0.28,0.67),log="y",main="",xlab='GC content',ylab="Read counts",bty ="n",yaxt='n')
axis(2,at=c(0,2^c(0,2,4,6,8))+1,labels=c(0,2^c(0,2,4,6,8)))


### gc variation ~ region width
library(BSgenome.Hsapiens.UCSC.hg19)
gcct <- GRangesList()
binsizes <- c(200L,300L,400L,500L)
for(i in seq_along(binsizes)){
    binsize <- binsizes[i]
    seqs <- seqlengths(Hsapiens)[1:24]
    starts <- lapply(seqs, function(i) seq(1, i, binsize))
    ends <- lapply(seqs, function(i) {
        if(i%%binsize == 0) seq(binsize, i, binsize)
        else c(seq(binsize, i, binsize), i)
    })
    chrs <- rep(names(seqs), times=sapply(starts, length))
    regions <- GRanges(chrs, IRanges(start=unlist(starts),
                                     end=unlist(ends)))
    seqlengths(regions) <- seqs[seqlevels(regions)]
    seqs <- getSeq(Hsapiens,regions)
    gcpos <- alphabetFrequency(seqs)
    gc <- rowSums(gcpos[,c('C','G')])/binsize
    mcols(regions)$gc <- gc
    gcct[[i]] <- regions
    save(gcct,file='~/share/gcapc/20160912.human.gcct.rda')
    cat(i,'\n')
}
for(i in seq_along(gcct)){
    mcols(gcct[[i]])$gc[mcols(gcct[[i]])$gc==0] <- NA
}
rgidx <- seq_along(gcct[[length(gcct)]])
gcdat <- matrix(0,length(rgidx),length(gcct))
for(i in seq_len(length(gcct)-1)){
    fo <- findOverlaps(gcct[[i]],gcct[[length(gcct)]])
    gcdat[,i] <- mcols(gcct[[i]])$gc[queryHits(fo)[match(rgidx,subjectHits(fo))]]
}
gcdat[,length(gcct)] <- mcols(gcct[[length(gcct)]])$gc
save(gcdat,file='~/share/gcapc/20160912.human.gcdat.rda')
par(cex.axis=1.2,cex.lab=1.2,cex.main=1.2,cex=1.1,lwd=2)
boxplot(gcdat,ylim=c(0.1,0.7),outline=F,names=binsizes,
        ylab="GC content",xlab="Genome region size (bp)",main="GC Content Distribution")


### gc variation ~ flanking width ## need information from previous plot
library(BSgenome.Hsapiens.UCSC.hg19)
bdwidth <- 200
seqs <- seqlengths(Hsapiens)[1:24]
seqs <- floor(seqs/bdwidth-2)*bdwidth
starts <- lapply(seqs, function(i) seq(1+bdwidth*2, i, bdwidth))
ends <- lapply(seqs, function(i) seq(bdwidth*3, i, bdwidth))
chrs <- rep(names(seqs), times=sapply(starts, length))
samp <- 0.05
sampidx <- sort(sample.int(length(chrs),ceiling(length(chrs)*samp)))
region <- GRanges(chrs[sampidx], IRanges(start=unlist(starts)[sampidx],
                                         end=unlist(ends)[sampidx]))
#region <- GRanges(chrs, IRanges(start=unlist(starts),end=unlist(ends)))
gccont0 <- function(region,flank=200L,genome="BSgenome.Hsapiens.UCSC.hg19"){
    library(GenomicRanges)
    library(genome,character.only=TRUE)
    cat("......... flanking:",flank,'\n')
    rwidth <- width(region[1])
    nr <- shift(resize(region,rwidth + flank*2),-flank)
    seqs <- getSeq(Hsapiens,nr)
    gcpos <- gregexpr('(G)|(C)',seqs)
    weight <- c(seq_len(flank),rep(flank+1,rwidth),rev(seq_len(flank)))
    gcnumw <- sapply(gcpos,function(x) sum(weight[x]))
    gcnumw[gcnumw==(sum(weight)-1)] <- 0
    round(gcnumw / (flank+rwidth) / (flank+1),3)
}
flanks <- c(0,50,100,150)
gcf <- matrix(0,length(region),length(flanks))
for(i in seq_along(flanks)){
    gcf[,i] <- gccont0(region,flank=flanks[i])
    cat(i,'\n')
    save(gcf,file="~/share/gcapc/20160913.human.gcf.region200.flank0-150.50.rda")
}
gcf[gcf==0] <- NA
layout(matrix(1:2,1,2))
s1 <- sample(seq_len(nrow(gcdat)),100000)
s2 <- sample(seq_len(nrow(gcf)),100000)
par(cex.axis=1.8,cex.lab=1.8,cex.main=1.8,cex=1.1)
plot(gcdat[s1,2],gcdat[s1,3],col=rgb(0,0,0,alpha=0.2),pch=20,ylim=c(0.2,0.8),xlim=c(0.2,0.8),
     xlab='',ylab="",main="",yaxt='n',xaxt='n')
plot(gcf[s2,2],gcf[s2,3],col=rgb(0,0,0,alpha=0.2),pch=20,ylim=c(0.2,0.8),xlim=c(0.2,0.8),
          xlab='',ylab="",main="",yaxt='n',xaxt='n')
layout(matrix(1:2,1,2))
par(cex.axis=1.2,cex.lab=1.2,cex.main=1.2,cex=1.1)
plot(NULL,xlab='300bp fragments',ylab="400bp fragments",main="GC Content",ylim=c(0.2,0.8),xlim=c(0.2,0.8),bty='n')
plot(NULL,NULL,xlab='300bp region (200bp in + 50*2 out)',ylab="400bp region (200bp in + 100*2 out)",main="Effective GC Content",
     ylim=c(0.2,0.8),xlim=c(0.2,0.8),bty='n')

