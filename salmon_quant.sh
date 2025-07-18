#!/bin/bash
# This file utilizes Salmon software (unix OSs only) to quantify transcript counts
# listed in fastq files in your working directory
# note that you can run Salmon on Windows but
# you must use a Linux virtual machine
# Salmon website: https://combine-lab.github.io/salmon/
# How to run Salmon on windows: https://www.protocols.io/view/ubuntu-on-windows-for-computational-biology-sfuebnw.pdf

files1=(data/*R1.fastq.gz) # all fastq files ending in R1 aka the first paired reads on that sample
files2=(data/*R2.fastq.gz) # all fastq files ending in R2 aka the second paired end reads
for i in {0..74}
do
	samp1=`basename -s .fastq.gz ${files1[i]}`
	samp2=`basename -s .fastq.gz ${files2[i]}`
	printf 'Processing samples %s and %s\n' ${samp1} ${samp2}
	salmon quant -i mm_index -l A --gcBias\
		-1 data/${samp1}.fastq.gz \
		-2 data/${samp2}.fastq.gz \
		-p 16 --validateMappings -o quants/${samp1}_quant

done