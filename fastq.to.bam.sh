#!/bin/bash
module load centos6/bwa-0.7.7
module load centos6/samtools-0.1.19

###### sequencing reads alignment using bwa
fastq=(*.fastq.*gz)
echo ${#fastq[@]}

for ((i=0; i<=${#fastq[@]}-1; i++))
do
    f=${fastq[$i]}
    echo $f
    n=(${f//./ })
    if [ ! -e ${n[0]}.bam ]; then
	echo aligning...
	bwa aln female.hg19.fa $f | \
	    bwa samse female.hg19.fa - $f | \
	    samtools view -Su - | \
	    samtools sort - ${n[0]}
    	samtools index ${n[0]}.bam
    fi
done