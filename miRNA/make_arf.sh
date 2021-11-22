#!/bin/sh

bam=/Users/ryan/Documents/GitHub/sm_RNA_seq/SRR1187947/Mapped_Reads/SRR1187947_trimmed.fq.gz.mapped.md.bam

conda activate mirdeep

samtools view -h $bam > SRR1187947_trimmed.fq.gz.mapped.md.sam

sam=SRR1187947_trimmed.fq.gz.mapped.md.sam


#Throws an error that there is not actually and output file

#bwa_sam_converter.pl -i $sam -a SRR1187947_trimmed.fq.gz.mapped.md.reads.arf

bwa_sam_converter.pl -i $sam -c -o reads_collapsed.fa -a reads_collapsed_vs_genome.arf