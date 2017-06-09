#!/bin/bash
module load centos6/MACS-2.0.10_python-2.7.3

###### call peaks using MACS2
bams=(*.bam)
folder=macs2
for ((i=0; i<=${#bams[@]}-1; i++))
do
    bam=${bams[i]}
    n=(${bam//./ })
    echo $n
    if [ ! -e $folder/$n\_peaks.xls ]; then
	macs2 callpeak -t $bam -f BAM -g hs -n $n -q 0.99 --outdir $folder
    fi
done

