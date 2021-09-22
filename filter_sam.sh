#!/bin/sh

#Filter sam for different alingments 


#View Uniqely Mapping Reads 
#samtools view -bq 1 file.bam > unique.bam

#Get only secondary mapped reads
samtools view -f 256 SRR1187947_trimmed.fq.gz.bowtie.multimapping_test.sam  > multi_mapped.sam