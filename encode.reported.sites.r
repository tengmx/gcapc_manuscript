library(GenomicRanges)
MARKER <- "CTCF"
marker <- list(CTCF="Ctcf")
celllst <- list(CTCF=c("Gm12878","Helas3","Hepg2","Huvec","K562","Nhek"))



###### site read count
library(Rsubread)
library(Rsamtools)
library(BSgenome.Hsapiens.UCSC.hg19)
### total sites from encode bed files, including proximal and distal
chipseqp <- read.table(paste0("data/track/chipseq_track_proximal.bed"),
                       sep='\t',header=F,stringsAsFactors=F)
chipseqd <- read.table(paste0("data/track/chipseq_track_distal.bed"),
                       sep='\t',header=F,stringsAsFactors=F)
sitep <- GRanges(chipseqp[,1],IRanges(start=chipseqp[,2]+1,end=chipseqp[,3]),
                   label=chipseqp[,4],detail=chipseqp[,6],tag="proximal")
sited <- GRanges(chipseqd[,1],IRanges(start=chipseqd[,2]+1,end=chipseqd[,3]),
                   label=chipseqd[,4],detail=chipseqd[,6],tag="distal")
sites <- sort(sortSeqlevels(c(sitep,sited)))
seqlevels(sites,force=T) <- intersect(seqlevels(sites),seqlevels(Hsapiens)[1:22])
seqlengths(sites) <- seqlengths(Hsapiens)[seqlevels(sites)]
### sites particularly for MARKER, e.g CTCF, POLR2A
markidx <- sapply(strsplit(mcols(sites)$detail,',|\\('),function(x) any(x %in% MARKER))
sitesDF <- as.data.frame(sites)
sitesDF <- data.frame(name=seq_len(nrow(sitesDF)),sitesDF[,c(1:3,5:6,8)])
colnames(sitesDF) <- c('GeneID', 'Chr', 'Start', 'End','Strand','Label','Tag')
### cell lines have data from all labs: in total 40 bam files for CTCF
bamfiles <- list.files(paste0("data/",marker[[MARKER]],"/all"),pattern="^wgEncode.+.bam$",full=T)
load(paste0("data/",marker[[MARKER]],"/all/meta.rda"))
cells <- celllst[[MARKER]]
meta <- meta[meta[,3] %in% cells,]
bamfiles <- bamfiles[match(gsub(".fastq.gz","",meta[,1]),gsub(".bam","",basename(bamfiles)))]
#identical(gsub(".bam","",basename(bamfiles)),gsub(".fastq.gz","",meta[,1]))
### sites raw reads count with each read extended to 400bp
site_raw <- featureCounts(bamfiles,annot.ext=sitesDF,minMQS=10,readExtension3=400-150-36,
                          allowMultiOverlap=TRUE,ignoreDup=TRUE)
colnames(site_raw$counts) <- gsub(paste0("(data.",marker[[MARKER]],".all.wgEncode)|(.bam)"),"",
                                           colnames(site_raw$counts))
colnames(site_raw$stat) <- gsub(paste0("(data.",marker[[MARKER]],".all.wgEncode)|(.bam)"),"",
                                           colnames(site_raw$stat))
names(sites) <- rownames(site_raw$counts)
#hist(log10(colSums(site_raw$counts[markidx,])),nclass=40) ## filter: sum>=1e+6
#filteridx <- which(colSums(site_raw$counts[markidx,]) >= 1e+6) #CTCF
filteridx <- which(colSums(site_raw$counts[markidx,]) >= 4e+5) #POLR2A
rawsig <- site_raw$counts[,filteridx]
rawstat <- site_raw$stat[,c(1,filteridx+1)]
meta <- meta[filteridx,]
save(rawsig,rawstat,sites,meta,markidx,
     file=paste0("rdas/",marker[[MARKER]],
         "/v3_sitereadcount/20160314.sites.raw38.nodup.rda"))


###### site GC, assume fragsize 400bp ## weighted GC, where the site is weighted more
load(paste0("rdas/",marker[[MARKER]],"/v3_sitereadcount/20160314.sites.raw38.nodup.rda")) 
library(BSgenome.Hsapiens.UCSC.hg19)
fragsize <- 400
sitegc <- c()
num <- 0
for(chr in seqlevels(sites)){
    startchr <- start(sites[seqnames(sites)==chr])
    endchr <- end(sites[seqnames(sites)==chr])
    num <- num + length(startchr)
    start <- lapply(seq_len(length(startchr)),function(i) (startchr[i] - fragsize + 150):(endchr[i] - 149))
    startvec <- unlist(start)
    view <- Views(Hsapiens[[chr]],start=startvec,width=fragsize)
    wsf <- alphabetFrequency(view, baseOnly=T, as.prob=T)
    gcf <- rowSums(wsf[,c("C","G")])/rowSums(wsf)
    gcf[wsf[,'other'] > 0.1] <- NA
    gcv <- Views(gcf,end=cumsum(sapply(start,length)),width=sapply(start,length))
    gc <- mean(gcv)
    sitegc <- c(sitegc,gc)
    cat(chr,length(startchr),num,'\n')
}
#identical(order(seqnames(sites)),seq_len(length(sites)))
save(sitegc,sites,fragsize,file="rdas/20160314.hg19.gc.150bpsites.400fragsize.rda")


###### mixture E-M model
ll <- function(y,mu1,mu0,p,z){
    logp1 <- dpois(y, lambda = mu1, log = TRUE)
    logp0 <- dpois(y, lambda = mu0, log = TRUE)
    sum(z*(logp1+log(p)) + (1-z)*(logp0+log(1-p)))
}
estep <- function(y,mu1,mu0,p){
    logp1 <- dpois(y, lambda = mu1, log = TRUE)
    logp0 <- dpois(y, lambda = mu0, log = TRUE)
    z <- 1/(1+exp(logp0-logp1)*(1-p)/p)  #(p1*p)/(p1*p+p0*(1-p))
    z
}
mstep1 <- function(y,z,gc){ # f(gc) + b
    library(splines)
    p <- (2+sum(z))/(2*2+length(z))
    dat <- data.frame(y=y,gc=gc)
    dat1 <- dat[z>=0.5,]
    dat0 <- dat[z<0.5,]
    lmns0 <- glm(y ~ ns(gc, df = 2), data=dat0, family="poisson")
    pred1Y <- predict(lmns0,data.frame(gc = dat1$gc),type="response")
    expB <- sum(dat1$y * pred1Y)/sum(pred1Y^2)
    predY <- predict(lmns0, data.frame(gc = gc),type="response")
    mu1 <- predY * expB
    mu0 <- predY
    list(p=p,mu1=mu1,mu0=mu0,z=z)
}
mstep <- mstep2 <- function(y,z,gc){ # f1(gc); f2(gc)
    library(splines)
    p <- (2+sum(z))/(2*2+length(z))
    dat <- data.frame(y=y,gc=gc)
    dat1 <- dat[z>=0.5,]
    dat0 <- dat[z<0.5,]
    lmns0 <- glm(y ~ ns(gc, df = 2), data=dat0, family="poisson")
    lmns1 <- glm(y ~ ns(gc, df = 2), data=dat1, family="poisson") #
    predY0 <- predict(lmns0, data.frame(gc = gc),type="response")
    predY1 <- predict(lmns1, data.frame(gc = gc),type="response") #
    list(p=p,mu1=predY1,mu0=predY0,z=z)
}
em <- function(y,gc,mu1=50,mu0=5,p=0.01){
    z <- estep(y,mu1,mu0,p)
    i <- 0
    llf <- ll(y,mu1,mu0,p,z)
    llgap <- llf
    while(abs(llgap) > (abs(llf) * 1e-5) && i < 100){
        para <- mstep(y,z,gc)
        z <- estep(y,para$mu1,para$mu0,para$p)
        if(sum(z>=0.5) < length(gc)*0.001 | sum(z<0.5) < length(gc)*0.001) break;
        lli <- ll(y,para$mu1,para$mu0,para$p,z)
        llgap <- lli - llf
        llf <- lli
        i <- i + 1
        cat(i,'\t',llf,'\t',llgap,'\n')
    }
    list(para=para,iter=i)
}
load(paste0("rdas/",marker[[MARKER]],"/v3_sitereadcount/20160314.sites.raw38.nodup.rda"))
load(paste0("rdas/20160314.hg19.gc.150bpsites.400fragsize.rda"))
layout(matrix(seq_len(ceiling(ncol(rawsig)/4)*4),ceiling(ncol(rawsig)/4),4,byrow=T))
par(cex.lab=2,cex.axis=2,cex.main=2,font.axis=2,font.lab=2) 
fitarg <- list()
for(rep in seq_len(ncol(rawsig))){
    args <- c(128,4,0.06,0.3,0.8,0)
    y0 <- rawsig[markidx,rep]
    gc0 <- sitegc[markidx]
    idx <- gc0>=args[4] & gc0<=args[5] & y0>=args[6] & !is.na(y0) & !is.na(gc0)
    y <- y0[idx]
    gc <- gc0[idx]
    ems <- em(y,gc,mu1=args[1],mu0=args[2],p=args[3])
    p <- ems$para$p
    mu1 <- mu0 <- z <- array(NA,dim=length(idx))
    mu1[idx] <- ems$para$mu1
    mu0[idx] <- ems$para$mu0
    z[idx] <- ems$para$z
    fitarg[[rep]] <- list(mu1=mu1,mu0=mu0,z=z,idx=idx,p=p)
    para <- ems$para
    iter <- ems$iter
    plot(gc,y+0.5,col=rgb(0,0,0,alpha=0.1),xlim=args[4:5],pch=20,
         main=paste0("rep: ",rep,"; iter: ",iter),xlab='GC',ylab="RC",log='y',yaxt='n')
    idx <- sample(seq_along(gc),1000)
    idx <- idx[order(gc[idx])]
    lines(gc[idx],para$mu1[idx]+0.5,col='red',lwd=3)
    lines(gc[idx],para$mu0[idx]+0.5,col='blue',lwd=3)
    axis(side=2, at=c(0,2^(0:10))+0.5, labels=c(0,2^(0:10)))
    cat(rep,"\tp:",para$p,"\n")
}
dev.off()
save(fitarg,file=paste0("rdas/",marker[[MARKER]],
                "/v3_sitereadcount/20160329.150bpsites.sigAndGc.mixture.fgcdiff.poisson.smallset.rda"))

pcaplot <- function(pcatbl,colors,pchs,main,legend=TRUE,first3=FALSE){
    library(corpcor)
    svd_tbl <- wt.scale(t(log2(pcatbl+1)),center=TRUE,scale=TRUE)
    svd_run <- fast.svd(svd_tbl)
    pcas <- svd_run$u
    props <- round(svd_run$d^2/sum(svd_run$d^2),3)
    library(RColorBrewer)
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
        plot(pcas[,1],pcas[,3],col=color,main=main,pch=pch,
             xlab=paste("PC1:",props[1]*100,"%"),ylab=paste("PC3:",props[3]*100,"%"))
        plot(pcas[,2],pcas[,3],col=color,main=main,pch=pch,
             xlab=paste("PC2:",props[2]*100,"%"),ylab=paste("PC3:",props[3]*100,"%"))
    }
}





###### plots
# 1a pca
###### distinguish cell type better
cells <- c("HUVEC","HeLa-S3","HepG2","GM12878","K562","NHEK")
load(paste0("rdas/",marker[[MARKER]],"/v3_sitereadcount/20160314.sites.raw38.nodup.rda"))
load(paste0("rdas/20160314.hg19.gc.150bpsites.400fragsize.rda"))
load(paste0("rdas/",marker[[MARKER]],"/v3_sitereadcount/20160329.150bpsites.sigAndGc.mixture.fgcdiff.poisson.smallset.rda"))
cellmarkidx <- grepl(paste0(MARKER,"\\(",cells[1],"\\)"), mcols(sites)$detail) |
    grepl(paste0(MARKER,"\\(",cells[2],"\\)"), mcols(sites)$detail) |
    grepl(paste0(MARKER,"\\(",cells[3],"\\)"), mcols(sites)$detail) |
    grepl(paste0(MARKER,"\\(",cells[4],"\\)"), mcols(sites)$detail) |
    grepl(paste0(MARKER,"\\(",cells[5],"\\)"), mcols(sites)$detail) |
    grepl(paste0(MARKER,"\\(",cells[6],"\\)"), mcols(sites)$detail)
sig <- rawsig[cellmarkidx,]
gc <- sitegc[cellmarkidx]
met <- meta
sit <- sites[cellmarkidx]
fit <- fitarg
for(i in seq_len(length(fit))){
    fit[[i]]$mu1 <- fit[[i]]$mu1[cellmarkidx[markidx]]
    fit[[i]]$mu0 <- fit[[i]]$mu0[cellmarkidx[markidx]]
    fit[[i]]$z <- fit[[i]]$z[cellmarkidx[markidx]]
    fit[[i]]$idx <- fit[[i]]$idx[cellmarkidx[markidx]]
}
sizefac <- colSums(sig)
sigbyscale <- t(t(sig)/sizefac*median(sizefac))
sigbygc <- sig
for(i in seq_len(length(fit))){
    gce <- log2((fit[[i]]$mu1/median(fit[[i]]$mu1[fit[[i]]$z>=0.5],na.rm=T))) * fit[[i]]$z +
           log2((fit[[i]]$mu0/median(fit[[i]]$mu0[fit[[i]]$z>=0.5],na.rm=T))) * (1-fit[[i]]$z)
    sigbygc[,i] <- sig[,i] / 2^gce
}
sizefac <- colSums(sigbygc,na.rm=T)
sigbygc <- t(t(sigbygc)/sizefac*median(sizefac))

pdf("~/share/gcapc/20161122.figure5.pdf",15,6)
layout(matrix(c(1,1,2,2,3),1,5))
par(cex.axis=1.2,cex.lab=1.2,font.lab=2,cex.main=1.2,lwd=2,cex=1.05)
idx <- !is.na(sigbygc[,1])
pcaplot(sigbyscale[idx,],meta[,3],meta[,4],main="With GC Effects",first3=FALSE)
pcaplot(sigbygc[idx,],meta[,3],meta[,4],main="Without GC Effects",first3=FALSE)

####################################################
cell1 <- "HUVEC" #c("HUVEC","HeLa-S3","HepG2","GM12878","K562","NHEK")
cell2 <- "Huvec" #c("Huvec","Helas3","Hepg2","Gm12878","K562","Nhek")
cellmarkidx <- grepl(paste0(MARKER,"\\(",cell1,"\\)"), mcols(sites)$detail)
sig <- rawsig[cellmarkidx,meta[,3]==cell2]
gc <- sitegc[cellmarkidx]
met <- meta[meta[,3]==cell2,]
sit <- sites[cellmarkidx]
fit <- fitarg[meta[,3]==cell2]
for(i in seq_len(length(fit))){
    fit[[i]]$mu1 <- fit[[i]]$mu1[cellmarkidx[markidx]]
    fit[[i]]$mu0 <- fit[[i]]$mu0[cellmarkidx[markidx]]
    fit[[i]]$z <- fit[[i]]$z[cellmarkidx[markidx]]
    fit[[i]]$idx <- fit[[i]]$idx[cellmarkidx[markidx]]
}
labfac <- as.factor(met[,4])
sizefac <- colSums(sig)
sigbyscale <- t(t(sig)/sizefac*median(sizefac))
sigbygc <- sig
for(i in seq_len(length(fit))){
    gce <- log2((fit[[i]]$mu1/median(fit[[i]]$mu1[fit[[i]]$z>=0.5],na.rm=T))) * fit[[i]]$z +
           log2((fit[[i]]$mu0/median(fit[[i]]$mu0[fit[[i]]$z>=0.5],na.rm=T))) * (1-fit[[i]]$z)
    sigbygc[,i] <- sig[,i] / 2^gce
}
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
save(meansqbyscale,meansqbygc,
              file=paste0("~/workspace/chipseq/gcapc/rdas/20161122.sites.",cell2,".scale.gccorrect.meansq.rda"))
load("~/workspace/chipseq/gcapc/rdas/20161122.sites.huvec.scale.gccorrect.meansq.rda")
boxplot(meansqbyscale,meansqbygc,main=cell1,ylim=c(0,5),outline=F,ylab="Mean square",
        names=c('With GC effects','Without GC effects'))

