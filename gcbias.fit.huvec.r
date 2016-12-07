library(gcapc)
library(BSgenome)
library(splines)
bams <- c('wgEncodeUwTfbsHuvecCtcfStdRawDataRep1.bam','wgEncodeUwTfbsHuvecCtcfStdRawDataRep2.bam',
          'wgEncodeOpenChromChipHuvecCtcfRawDataRep1.bam','wgEncodeOpenChromChipHuvecCtcfRawDataRep2.bam',
          'wgEncodeBroadHistoneHuvecCtcfStdRawDataRep1.bam','wgEncodeBroadHistoneHuvecCtcfStdRawDataRep2.bam')
for(bam in seq_along(bams)){
    cov <- read5endCoverage(bams[bam],chroms=paste0('chr',1:22))
    bdwidth <- bindWidth(cov, range = c(50L, 400L), step = 50L, odd = TRUE)
    samp = 0.05;gcrange = c(0.3, 0.8); mu0 = 1; mu1 = 50; p = 0.02; converge = 0.001
    genome = "hg19";emtrace=TRUE
    genome <- getBSgenome(genome)
    bdw <- bdwidth[1]
    halfbdw <- floor(bdw/2)
    pdwh <- bdwidth[2]
    flank <- pdwh-bdw+halfbdw
    ### regions and reads count
    cat("Starting to estimate GC effects.\n")
    cat("...... Sampling regions\n")
    seqs <- sapply(cov$fwd,length)
    seqs <- floor(seqs/bdw-2)*bdw
    starts <- lapply(seqs, function(i) seq(1+bdw*2, i, bdw))
    ends <- lapply(seqs, function(i) seq(bdw*3, i, bdw))
    chrs <- rep(names(seqs), times=sapply(starts, length))
    sampidx <- sort(sample.int(length(chrs),ceiling(length(chrs)*samp)))
    region <- GRanges(chrs[sampidx], IRanges(start=unlist(starts)[sampidx],
                                             end=unlist(ends)[sampidx]))
    cat("......... Estimating using",length(region),"regions\n")
    regionsp <- resize(split(region,seqnames(region)),pdwh)
    cat("...... Counting reads\n")
    rcfwd <- unlist(viewSums(Views(cov$fwd,
                                   ranges(shift(regionsp,-flank)))))
    rcrev <- unlist(viewSums(Views(cov$rev,
                                   ranges(shift(regionsp,halfbdw)))))
    ### effective gc content
    cat("...... Calculating GC content with flanking",flank,"\n")
    rwidth <- width(region[1])
    nr <- shift(resize(region,rwidth + flank*2),-flank)
    seqs <- getSeq(genome,nr)
    gcpos <- startIndex(vmatchPattern("S", seqs, fixed="subject"))
    weight <- c(seq_len(flank),rep(flank+1,bdw),rev(seq_len(flank)))
    weight <- weight/sum(weight)
    gc <- round(sapply(gcpos,function(x) sum(weight[x])),3)
    ### em algorithms
    cat("...... Estimating GC effects\n")
    gc <- rep(gc,2)
    rc <- c(rcfwd,rcrev)
    idx <- gc>=gcrange[1] & gc<=gcrange[2] & !is.na(rc) & !is.na(gc)
    y <- rc[idx]
    gc <- gc[idx]
    logp1 <- dpois(y, lambda = mu1, log = TRUE)
    logp0 <- dpois(y, lambda = mu0, log = TRUE)
    z <- 1/(1+exp(logp0-logp1)*(1-p)/p)
    llf <- sum(z*(logp1+log(p)) + (1-z)*(logp0+log(1-p)))
    llgap <- llf
    i <- 0
    while(abs(llgap) > (abs(llf) * converge) && i < 100){
        p <- (2+sum(z))/(2*2+length(z))
        dat <- data.frame(y=y,gc=gc)
        dat1 <- dat[z>=0.5,]
        dat0 <- dat[z<0.5,]
        lmns0 <- glm(y ~ ns(gc, df = 2), data=dat0, family="poisson")
        lmns1 <- glm(y ~ ns(gc, df = 2), data=dat1, family="poisson")
        predY0 <- predict(lmns0, data.frame(gc = gc),type="response")
        predY1 <- predict(lmns1, data.frame(gc = gc),type="response")
        logp1 <- dpois(y, lambda = predY1, log = TRUE)
        logp0 <- dpois(y, lambda = predY0, log = TRUE)
        z <- 1/(1+exp(logp0-logp1)*(1-p)/p)
        if(sum(z>=0.5) < length(gc)*0.0005 | sum(z<0.5) < length(gc)*0.0005)
            break;
        lli <- sum(z*(logp1+log(p)) + (1-z)*(logp0+log(1-p)))
        llgap <- lli - llf
        llf <- lli
        i <- i + 1
        if(emtrace)
            cat("......... Iteration",i,'\tll',llf,'\tincrement',llgap,'\n')
    }
    ### gc effects
    gcbias <- list(glm0=lmns0,glm1=lmns1,
                   mu0=predY0,mu1=predY1,z=z,p=p,ll=llf,y=y,gc=gc,gcrange=gcrange,bdwidth=bdwidth)
    save(gcbias,file=paste0('~/workspace/chipseq/gcapc/rdas/20161122.',bams[bam],'.gcbias.rda'))
}

library(RColorBrewer)
layout(matrix(1:6,2,3))
par(cex.axis=1.8,cex.lab=1.8,cex.main=1.8,cex=1.1,lwd=2)
rcb <- rev(brewer.pal(5,"RdYlBu"))
rbPal <- colorRampPalette(rcb)
for(bam in bams){
    load(paste0('~/workspace/chipseq/gcapc/rdas/20161122.',bam,'.gcbias.rda'))
    tmp <- length(gcbias$y)
    idx0 <- sample.int(tmp,min(150000,tmp))
    color <- rbPal(5)[as.numeric(cut(gcbias$z[idx0],breaks = 5))]
    plot(gcbias$gc[idx0],gcbias$y[idx0]+0.5,col=color,xlim=gcbias$gcrange,ylim=c(0.5,80),
         pch=20,main="",xlab='',ylab="",log='y',yaxt='n',xaxt='n')
    idx00 <- sample.int(tmp,min(5000,tmp))
    idx00 <- idx00[order(gcbias$gc[idx00])]
    lines(gcbias$gc[idx00],gcbias$mu1[idx00]+0.5,col=rcb[5],lwd=10,lty=2)
    lines(gcbias$gc[idx00],gcbias$mu0[idx00]+0.5,col=rcb[1],lwd=10,lty=2)
}
par(cex.axis=1.2,cex.lab=1.2,cex.main=1.2,cex=1.1,lwd=2)
plot(NULL,ylim=c(0.5,80),xlim=gcbias$gcrange,log="y",main="",
     xlab='Effective GC content',ylab="Read counts",yaxt='n',bty ="n")
axis(side=2, at=c(0,2^(0:10))+0.5, labels=c(0,2^(0:10)))
legend('topleft',legend=c('0-0.2','0.2-0.4','0.4-0.6','0.6-0.8','0.8-1'),pch=20,col=rcb,bty='n')

