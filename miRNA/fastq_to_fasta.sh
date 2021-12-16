#!/bin/sh

fastq=/Users/ryan/Documents/GitHub/sm_RNA_seq/SRR1187947/Trimmed_Reads/SRR1187947_trimmed.fq

#Convert the Raw Fastq File to a FASTA file (if needed)

gunzip -c $fastq > SRR1187947_trimmed.fq 

#Remove the quality scores
awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' SRR1187947_trimmed.fq > SRR1187947_trimmed.fasta

#remove the white space from the identifier

perl -i -pe 's|[ \t]||g' SRR1187947_trimmed.fasta
