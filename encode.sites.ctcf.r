###### analysis baed on CTCF ENCODE sites

library(GenomicRanges)
library(Rsubread)
library(BSgenome.Hsapiens.UCSC.hg19)
library(gcapc) ## version 1.0.8
library(corpcor)
library(RColorBrewer)
library(sva)

## load bamfiles and corresponding meta information: lab, cell, replicate
MARKER <- "CTCF"
marker <- "Ctcf"
CELLS <- c("HUVEC","HeLa-S3","HepG2","GM12878","K562","NHEK")
cells <- c("Huvec","Helas3","Hepg2","Gm12878","K562","Nhek")
bamfiles <- list.files(".",pattern="^wgEncode.+.bam$",full=T)
load(paste0("meta.rda"))
meta <- meta[meta[,3] %in% cells,]
bamfiles <- bamfiles[match(gsub(".fastq.gz","",meta[,1]),gsub(".bam","",basename(bamfiles)))]

## reads count by assuming fragment length as 400bp (equivalent to flank as 250bp + sites 150)
chipseqp <- read.table(paste0("data/track/chipseq_track_proximal.bed"),
                       sep='\t',header=F,stringsAsFactors=F)
chipseqd <- read.table(paste0("data/track/chipseq_track_distal.bed"),
                       sep='\t',header=F,stringsAsFactors=F)
sitep <- GRanges(chipseqp[,1],IRanges(start=chipseqp[,2]+1,end=chipseqp[,3]),
                   label=chipseqp[,4],detail=chipseqp[,6],tag="proximal")
sited <- GRanges(chipseqd[,1],IRanges(start=chipseqd[,2]+1,end=chipseqd[,3]),
                   label=chipseqd[,4],detail=chipseqd[,6],tag="distal")
sites <- sort(sortSeqlevels(c(sitep,sited)))
seqlevels(sites,force=T) <- intersect(seqlevels(sites),seqlevels(Hsapiens)[1:23])
seqlengths(sites) <- seqlengths(Hsapiens)[seqlevels(sites)]
markidx <- sapply(strsplit(mcols(sites)$detail,',|\\('),function(x) any(x %in% MARKER))
sites <- sites[markidx]
sitesDF <- as.data.frame(sites)
sitesDF <- data.frame(name=seq_len(nrow(sitesDF)),sitesDF[,c(1:3,5:6,8)])
colnames(sitesDF) <- c('GeneID', 'Chr', 'Start', 'End','Strand','Label','Tag')
siteraw <- featureCounts(bamfiles,annot.ext=sitesDF,minMQS=30,readExtension3=400-150-36,
                         minOverlap=1,primaryOnly=TRUE,allowMultiOverlap=TRUE,ignoreDup=TRUE)
colnames(siteraw$counts) <- gsub(paste0("(data.",marker,".all.wgEncode)|(.bam)"),"",
                                           colnames(siteraw$counts))
colnames(siteraw$stat) <- gsub(paste0("(data.",marker,".all.wgEncode)|(.bam)"),"",
                                           colnames(siteraw$stat))
names(sites) <- rownames(siteraw$counts)
rawsig <- siteraw$counts
rawstat <- siteraw$stat

## refine site signals
cellmarkidx <- grepl(paste0(MARKER,"\\(",CELLS[1],"\\)"), mcols(sites)$detail) |
    grepl(paste0(MARKER,"\\(",CELLS[2],"\\)"), mcols(sites)$detail) |
    grepl(paste0(MARKER,"\\(",CELLS[3],"\\)"), mcols(sites)$detail) |
    grepl(paste0(MARKER,"\\(",CELLS[4],"\\)"), mcols(sites)$detail) |
    grepl(paste0(MARKER,"\\(",CELLS[5],"\\)"), mcols(sites)$detail) |
    grepl(paste0(MARKER,"\\(",CELLS[6],"\\)"), mcols(sites)$detail)
newsig <- refineSites(rawsig,sites,flank=250L,outputidx=cellmarkidx,model="nbinom",converge=1e-3)

## PCA analysis
pcaplot <- function(pcatbl,colors,pchs,main,legend=TRUE,first3=FALSE,log=FALSE){
    if(!log){
        svd_tbl <- wt.scale(t(log2(pcatbl+1)),center=TRUE,scale=TRUE)
    }else{
        svd_tbl <- wt.scale(t(pcatbl),center=TRUE,scale=TRUE)
    }
    svd_run <- fast.svd(svd_tbl)
    pcas <- svd_run$u
    props <- round(svd_run$d^2/sum(svd_run$d^2),3)
    cells <- unique(colors)
    rcb <- brewer.pal(length(cells),'Dark2')
    names(rcb) <- cells
    pch <- as.numeric(as.factor(pchs))
    plot(pcas[,1],pcas[,2],col=rcb[colors],main=main,pch=pch,cex=1.2,
         xlab=paste("PC1:",props[1]*100,"%"),ylab=paste("PC2:",props[2]*100,"%"))
    if(legend){
        legend('bottomleft',cells,pch=1,col=rcb[cells],bty='n')
        legend('bottomright',unique(pchs),pch=unique(pch),bty='n')
    }
    if(first3){
        plot(pcas[,1],pcas[,3],col=rcb[colors],main=main,pch=pch,
             xlab=paste("PC1:",props[1]*100,"%"),ylab=paste("PC3:",props[3]*100,"%"))
        plot(pcas[,2],pcas[,3],col=rcb[colors],main=main,pch=pch,
             xlab=paste("PC2:",props[2]*100,"%"),ylab=paste("PC3:",props[3]*100,"%"))
    }
}
sig <- rawsig[cellmarkidx,]
sizefac <- colSums(sig,na.rm=T)
sigbyscale <- t(t(sig)/sizefac*median(sizefac))
sizefac <- colSums(newsig,na.rm=T)
sigbygc <- t(t(newsig)/sizefac*median(sizefac))
idx <- !is.na(sigbygc[,1])
pcaplot(sigbyscale[idx,],meta[,3],meta[,4],main="With GC Effects",first3=F)
pcaplot(sigbygc[idx,],meta[,3],meta[,4],main="Without GC Effects",first3=F)

## mean square
CELL <- "HUVEC"
cell <- "Huvec"
cell1markidx <- grepl(paste0(MARKER,"\\(",CELL,"\\)"), mcols(sites)$detail)
sig <- rawsig[cell1markidx,meta[,3]==cell]
sigbygc <- newsig[cell1markidx[cellmarkidx],meta[,3]==cell]
labfac <- as.factor(meta[meta[,3]==cell,4])
sizefac <- colSums(sig,na.rm=T)
sigbyscale <- t(t(sig)/sizefac*median(sizefac))
sizefac <- colSums(sigbygc,na.rm=T)
sigbygc <- t(t(sigbygc)/sizefac*median(sizefac))
meansq <- function(sig,fac){
    meansq <- array(NA,dim=nrow(sig))
    y <- log2(sig+0.5)
    for(i in seq_len(nrow(sig))){
        if(sum(is.na(y[i,]))==0)
            meansq[i] <- anova(lm(y[i,]~fac))$`Mean Sq`[1]
        cat(i,'\n')
    }
    meansq
}
meansqbyscale <- meansq(sigbyscale,labfac)
meansqbygc <- meansq(sigbygc,labfac)
boxplot(meansqbyscale,meansqbygc,main=CELL,ylim=c(0,5.3),outline=F,ylab="Cross laboratory variability",
        names=c('Before GC content adjustment','After GC content adjustment'))

## compare with combat
logsig <- log2(sig+1)
logsig <- logsig[rowSums(logsig==0)<ncol(sig),]
modcombat <- model.matrix(~1, data=data.frame(meta))
cbsig <- ComBat(dat=logsig, batch=as.factor(meta[,4]), mod=modcombat)
cbsig1 <- t(t(cbsig)-colMeans(cbsig)+median(colMeans(cbsig)))
idx <- !is.na(sigbygc[,1])
pcaplot(sigbygc[idx,],meta[,3],meta[,4],main="Our approach",first3=F)
pcaplot(cbsig1,meta[,3],meta[,4],main="After ComBat",first3=F,log=TRUE)
