###### run idr analysis for gcapc and SPP, then compare GC content and overlap
## IDR version 2.0.2

## peak files
library(GenomicRanges)
pfs <- list.files('.',pattern='^.+[1-2].peak$')
labs <- rep(c('Broad','UTA','UW'),each=2)

## run idr
for(i in 1:5){
  for(j in (i+1):6){
    if(labs[i]==labs[j]){
      system(paste0("~/anaconda/bin/idr --samples ",pfs[i]," ",pfs[j],
                    " --input-file-type bed --output-file sample",i,"_sample",j,"_",labs[i],"_",labs[j],".idr",
                    " --idr-threshold 0.02 --plot --rank 7 --peak-merge-method sum --initial-mu 0.2 --initial-sigma 1 ",
                    "--initial-rho 0.2 --initial-mix-param 0.5"))
      cat(i,j,'\n')
    }
  }
}

## idr peaks
idrfiles <- list.files('.',pattern='^sample.+idr$')
cellpeaklist <- GRangesList()
for(i in seq_along(idrfiles)){
  tmp <- read.table(idrfiles[i],sep='\t',header=F,stringsAsFactors=F)
  cellpeaklist[[i]] <- makeGRangesFromDataFrame(tmp,seqnames.field="V1",start.field="V2",end.field="V3",
                                                strand.field="V6",starts.in.df.are.0based=TRUE,keep.extra.columns=FALSE)
  cat(i, '\n')
}
lab <- labs[c(1,3,5)]
cellpeaks <- reduce(do.call(c,cellpeaklist))
fo <- lapply(cellpeaklist,function(x) findOverlaps(x,cellpeaks))
foidxlst <- lapply(fo,function(x) unique(subjectHits(x)))
foidx <- table(unlist(foidxlst))

## GC content
gccont.human <- function(grg){
  library(BSgenome.Hsapiens.UCSC.hg19)
  grgseq <- getSeq(Hsapiens,grg)
  wsf <- alphabetFrequency(grgseq, baseOnly=T, as.prob=T)
  rowSums(wsf[,c("C","G")])/rowSums(wsf)
}
gc <- gccont.human(cellpeaks)
gcs <- lapply(1:3,function(i) gc[intersect(which(foidx==1),foidxlst[[i]])])
boxplot(c(gcs,list(gc[foidx==3])),names=c(paste("Only",lab),"All labs"),
        ylab='GC content',pch=20,cex=0.3,ylim=c(0.2,0.86))

## venn diagram
library(VennDiagram)
gcs1 <- sapply(1:3,function(i) length(intersect(which(foidx==1),foidxlst[[i]])))
gcs2 <- sapply(1:3,function(i) length(intersect(intersect(which(foidx==2),foidxlst[[i]]),foidxlst[[i%%3+1]])))
gcs3 <- length(which(foidx==3))
grid.newpage()
draw.triple.venn(area1 = gcs3+gcs1[3]+gcs2[2]+gcs2[3],
                 area2 = gcs3+gcs1[2]+gcs2[1]+gcs2[2],
                 area3 = gcs3+gcs1[1]+gcs2[1]+gcs2[3], n12 = gcs3+gcs2[2], n23 = gcs3+gcs2[1], n13 = gcs3+gcs2[3],
                 n123 = gcs3, category = c("UW", "UTA", "Broad"), lwd=2,
                 fill = NULL)

## proportion of only by one lab 
sum(gcs1)/sum(gcs1,gcs2,gcs3)

