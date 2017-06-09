###### call peaks using gcapc, default parameters for CTCF

library(GenomicRanges)
library(gcapc) ## version 1.0.8
bams <- list.files(".",pattern='.+bam$')
for(reads in bams){
    sampling <- c(0.05,1) ## or c(0.15,1), c(0.25,1), c(0.05,5)
    mu0 <- 1 ## other TFs: 0.5
    mu1 <- 50 ## YY1 & EP300: 20; GATA2: 10
    theta0 <- 1 ## other TFs: 0.5
    theta1 <- 10 ## other TFs: 5
    prefilter <- 4 
    gcrange <- c(0.3,0.8) ## YY1 & EP300: c(0.25.0.75); GATA2: c(0.4,0.7)
    supervise <- FALSE ## or TRUE when supervised peak calling
    cov <- read5endCoverage(reads,chroms=paste0("chr",c(1:22,"X")))
    bdwidth <- bindWidth(cov,range = c(50L, 400L))
    if(supervise){
        tmp <- read.table(paste0(gsub(".bam","",reads),".peak"))
        supervise <- GRanges(tmp$V1,IRanges(tmp$V2+1,tmp$V3),sig=tmp$V7)
        gcbias <- gcEffects(cov,bdwidth,plot=T,sampling=sampling,model="nbinom",supervise=supervise,gcrange=gcrange,
                            mu0=mu0,mu1=mu1,theta0=theta0,theta1=theta1,converge=1e-3,gctype="ladder")
    }else{
        gcbias <- gcEffects(cov,bdwidth,plot=T,sampling=sampling,model="nbinom",gcrange=gcrange,
                            mu0=mu0,mu1=mu1,theta0=theta0,theta1=theta1,converge=1e-3,gctype="ladder")
    }
    peaks <- gcapcPeaks(cov,gcbias,bdwidth,plot=T,pv=0.99,gctype="ladder",prefilter=prefilter)
    save(reads,bdwidth,gcbias,peaks,mu0,mu1,theta0,theta1,sampling,prefilter,file=paste0(strsplit(reads,"\\.")[[1]][1],".rda"))
}

