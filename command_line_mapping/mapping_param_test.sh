#!/bin/sh 
#This script was used to check number of reads aligned and run time of bowtie with different presets

mkdir bowtie_performance_checks
cd bowtie_performance_checks
now=$(date)
echo $date > sensitive_start
bowtie2 --very-sensitive -a -p 4 -x /Volumes/BlytheLab_Files/HTSeq/Bowtie_Indices/dm6 -U /Users/ryan/Documents/GitHub/sm_RNA_seq/SRR1187947/Trimmed_Reads/SRR1187947_trimmed.fq.gz -S SRR1187947_mapped_verysensitive.sam
now=$(date) 
echo $date > sensitive_end 


#Note this code is slightly differen than above and is set to spit out a BAM This is for downstream process that I wanted to do with the file 
now=$(date) 
echo $date > sensitive_local_start 
bowtie2 --very-sensitive-local -a -p 4 -x /Volumes/BlytheLab_Files/HTSeq/Bowtie_Indices/dm6 -U /Users/ryan/Documents/GitHub/sm_RNA_seq/SRR1187947/Trimmed_Reads/SRR1187947_trimmed.fq.gz 2> SRR1187947_mapped_verysensitive_local.log.txt |
samtools view -bS - > RR1187947_mapped_verysensitive_local.mapped.bam
#-S SRR1187947_mapped_verysensitive_local.sam 

now=$(date)
echo $date > sensitive_local_end 

