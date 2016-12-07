library(gcapc)
library(BSgenome)
library(splines)
library(GenomicRanges)
bams <- c('wgEncodeUwTfbsHuvecCtcfStdRawDataRep1.bam','wgEncodeUwTfbsHuvecCtcfStdRawDataRep2.bam',
          'wgEncodeOpenChromChipHuvecCtcfRawDataRep1.bam','wgEncodeOpenChromChipHuvecCtcfRawDataRep2.bam',
          'wgEncodeBroadHistoneHuvecCtcfStdRawDataRep1.bam','wgEncodeBroadHistoneHuvecCtcfStdRawDataRep2.bam')
for(bam in seq_along(bams)){
    cov <- read5endCoverage(bams[bam],chroms=paste0('chr',1:22))
    bdwidth <- bindWidth(cov, range = c(50L, 400L), step = 50L, odd = TRUE)
    load(paste0('~/workspace/chipseq/gcapc/rdas/20161122.',bams[bam],'.gcbias.rda'))

    prefilter = 4L
    permute = 5L
    pv = 0.05
    genome <- getBSgenome('hg19')
    bdw <- bdwidth[1]
    halfbdw <- floor(bdw/2)
    pdwh <- bdwidth[2]
    flank <- pdwh - bdw + halfbdw
    cat("Starting to call peaks.\n")
    cat("...... prefiltering regions\n")
    seqs <- sapply(cov$fwd, length)
    seqs <- floor(seqs/bdw - 2) * bdw
    starts <- lapply(seqs, function(i) seq(1 + bdw * 2, i, bdw))
    ends <- lapply(seqs, function(i) seq(bdw * 3, i, bdw))
    chrs <- rep(names(seqs), times = sapply(starts, length))
    region <- GRanges(chrs, IRanges(start = unlist(starts), end = unlist(ends)))
    regionsp <- resize(split(region, seqnames(region)), pdwh)
    rcfwd <- unlist(viewSums(Views(cov$fwd, ranges(shift(regionsp,
                                                         -flank)))))
    rcrev <- unlist(viewSums(Views(cov$rev, ranges(shift(regionsp,
                                                         halfbdw)))))
    regions <- region[rcfwd + rcrev >= prefilter]
    rm(seqs, starts, ends, chrs, region, regionsp, rcfwd, rcrev)
    regionsrc <- reduce(shift(resize(regions, halfbdw * 6 + flank *
                                     2 + 1), -halfbdw * 2 - flank))
    regionsgc <- shift(resize(regionsrc, width(regionsrc) + halfbdw *
                              2), -halfbdw)
    rm(regions)
    cat("...... caculating GC content\n")
    nr <- shift(resize(regionsgc, width(regionsgc) + flank *
                       2), -flank)
    seqs <- getSeq(genome, nr)
    gcpos <- startIndex(vmatchPattern("S", seqs, fixed = "subject"))
    weight <- c(seq_len(flank), rep(flank + 1, bdw), rev(seq_len(flank)))
    weight <- weight/sum(weight)
    gcposb <- vector("integer", sum(width(nr)))
    gcposbi <- rep(seq_along(nr), times = width(nr))
    gcposbsp <- split(gcposb, gcposbi)
    for (i in seq_along(nr)) {
        gcposbsp[[i]][gcpos[[i]]] <- 1L
    }
    gcposbsprle <- as(gcposbsp, "RleList")
    gcnuml <- width(regionsgc) - halfbdw * 2
    gcnum <- vector("numeric", sum(gcnuml))
    gcnumi <- rep(seq_along(nr), times = gcnuml)
    gc <- split(gcnum, gcnumi)
    for (i in seq_along(nr)) {
        gc[[i]] <- round(runwtsum(gcposbsprle[[i]], k = length(weight),
                                  wt = weight), 3)
        if (i%%500 == 0)
            cat(".")
    }
    cat("\n")
    rm(regionsgc, nr, seqs, gcpos, gcposb, gcposbi, gcposbsp,
       gcposbsprle, gcnuml, gcnum, gcnumi)
    cat("...... caculating GC effect weights\n")
    gcbase <- round(seq(0, 1, 0.001), 3)
    mu0 <- predict(gcbias$glm0, data.frame(gc = gcbase), type = "response")
    gcwbase0 <- round(Rle(median(gcbias$mu0[gcbias$z > 0.5])/mu0),3)
    gcw <- RleList(lapply(gc, function(x) gcwbase0[x * 1000 +1]))
    rm(gc, gcbase, mu0, gcwbase0)
    regionrcsp <- split(regionsrc, seqnames(regionsrc))
    rcfwd <- RleList(unlist(viewApply(Views(cov$fwd, ranges(regionrcsp)),
                                      runsum, k = pdwh)))
    rcrev <- RleList(unlist(viewApply(Views(cov$rev, ranges(regionrcsp)),
                                      runsum, k = pdwh)))
    cat("...... estimating enrichment score\n")
    esl <- width(regionsrc) - halfbdw * 2 - flank * 2
    ir1 <- IRangesList(start = IntegerList(as.list(rep(1, length(esl)))),
                       end = IntegerList(as.list(esl)))
    ir2 <- IRangesList(start = IntegerList(as.list(rep(1 + halfbdw +
                           flank, length(esl)))), end = halfbdw + flank + IntegerList(as.list(esl)))
    ir3 <- IRangesList(start = IntegerList(as.list(rep(1 + halfbdw *
                           2 + flank * 2, length(esl)))), end = halfbdw * 2 + flank *
                       2 + IntegerList(as.list(esl)))
    rc1 <- rcfwd[ir1]
    rc2 <- rcfwd[ir2]
    rc3 <- rcrev[ir1]
    rc4 <- rcrev[ir2]
    gcw1 <- gcw[ir2]
    gcw2 <- gcw[ir3]
    gcw3 <- gcw[ir1]
    esrlt <- round(2 * sqrt(rc1 * rc4) * gcw1 - rc3 * gcw3 -
                   rc2 * gcw2, 3)
    names(esrlt) <- seq_along(regionsrc)
    rm(rcfwd, rcrev, gcw, rc1, rc2, rc3, rc4)
    cat("...... permutation analysis\n")
        esprlt <- RleList()
    for (i in seq_len(permute)) {
        covfwdp <- cov$fwd[ranges(regionrcsp)]
        covrevp <- cov$rev[ranges(regionrcsp)]
        for (i in seq_along(cov$fwd)) {
            covfwdp[[i]] <- covfwdp[[i]][sample.int(length(covfwdp[[i]]))]
            covrevp[[i]] <- covrevp[[i]][sample.int(length(covrevp[[i]]))]
        }
        end <- cumsum(width(regionrcsp))
        start <- end - width(regionrcsp) + 1
        regionrcspp <- IRangesList(start = start, end = end)
        rcfwdp <- RleList(unlist(viewApply(Views(covfwdp, regionrcspp),
                                           runsum, k = pdwh)))
        rcrevp <- RleList(unlist(viewApply(Views(covrevp, regionrcspp),
                                           runsum, k = pdwh)))
        rcp1 <- rcfwdp[ir1]
        rcp2 <- rcfwdp[ir2]
        rcp3 <- rcrevp[ir1]
        rcp4 <- rcrevp[ir2]
        esprlt <- c(esprlt, round(2 * sqrt(rcp1 * rcp4) * gcw1 -
                                  rcp3 * gcw3 - rcp2 * gcw2, 3))
    }
    perm <- as.numeric(unlist(esprlt, use.names = FALSE))
    rm(covfwdp, covrevp, start, end, regionrcspp, rcp1, rcp2,
       rcp3, rcp4, rcfwdp, rcrevp, gcw1, gcw2, gcw3, esprlt,
       esl, regionrcsp)
    cat("...... reporting peaks\n")
    sccut <- quantile(perm, 1 - pv)
    cat("......... enrichment scores cut at", sccut, "\n")
    save(esrlt,perm,pv,sccut,permute,prefilter,
         file=paste0('~/workspace/chipseq/gcapc/rdas/20161122.',bams[bam],'.perm.rda'))
    cat(bam,'\n')
}

layout(matrix(1:6,2,3))
par(cex.lab=1.2,cex.axis=1.2,lwd=2)
library(RColorBrewer)
colors <- brewer.pal(3,'RdGy')
for(bam in seq_along(bams)){
    load(paste0('~/workspace/chipseq/gcapc/rdas/20161122.',bams[bam],'.perm.rda'))
    es <- as.numeric(unlist(esrlt, use.names = FALSE))
    plot(density(perm, bw = 1), col = colors[3], xlab = "Enrichment score",
         xlim = c(-100,100), main = bams[bam])
    lines(density(es, bw = 1), col = colors[1])
    abline(v = sccut, lty = 2, col='grey')
    legend("topright", c("real", "perm"), lty = 1, col = colors[c(1,3)], bty = "n")
    cat(bam,'\n')
}

