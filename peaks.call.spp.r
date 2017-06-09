###### call peaks using SPP

library(spp) ## version 1.14
library(Rsamtools)
bams <- list.files(".",pattern='.+bam$')
fdr <- 0.99
for(bam in bams){
    ### read file
    chip.data <- read.bam.tags(bam)
    binding.characteristics <- get.binding.characteristics(chip.data,srange=c(80,500),bin=5)
    print(paste("binding peak separation distance =",binding.characteristics$peak$x))
    ### plot features
    par(mar = c(3.5,3.5,1.0,0.5), mgp = c(2,0.65,0), cex = 0.8)
    plot(binding.characteristics$cross.correlation,type='l',xlab="strand shift",ylab="cross-correlation")
    abline(v=binding.characteristics$peak$x,lty=2,col=2)
    ### prefilter tags
    chip.data <- select.informative.tags(chip.data,binding.characteristics)
    chip.data <- remove.local.tag.anomalies(chip.data)
    ### detect binding
    detection.window.halfsize <- binding.characteristics$whs
    bp <- find.binding.positions(signal.data=chip.data,fdr=fdr,whs=detection.window.halfsize)
    print(paste("detected",sum(unlist(lapply(bp$npl,function(d) length(d$x)))),"peaks"))
    ### output peaks
    write.narrowpeak.binding(bp,paste0("spp/",strsplit(bam,"\\.")[[1]][1],".peak"),
                             margin = round(detection.window.halfsize/2),npeaks=2000000)
    save(bp,binding.characteristics,file=paste0("spp/",strsplit(bam,"\\.")[[1]][1],".rda"))
}



