#!/bin/bash
module load centos6/bwa-0.7.7
module load centos6/samtools-0.1.19

## M: Hepg2
## F: K562, Gm12878, Helas3
## U: Huvec, Nhek

cell="Hepg2"
cd ~/workspace/chipseq/chipseqvar/encode2/data
cd $cell
if [ $cell = "Hepg2" ] 
then 
    sex="male" 
else
    sex="female" 
fi
fastq=$( ls *.fastq.gz )
for f in ${fastq[@]}
do
    echo $f
    n=(${f//./ })
    echo ${n[0]}
    echo $sex.hg19.fa
    bwa aln ~/g/hg19/encodedcc/referenceSequences/$sex.hg19.fa $f | bwa samse ~/g/hg19/encodedcc/referenceSequences/$sex.hg19.fa - $f | samtools view -Su - | samtools sort - ${n[0]}
done