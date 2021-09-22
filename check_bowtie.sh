#!/bin/sh 

mkdir bowtie_performance_checks
cd bowtie_performance_checks
now=$(date)
echo $date > sensitive_start
bowtie2 --very-sensitive -a -p 4 -x /Volumes/BlytheLab_Files/HTSeq/Bowtie_Indices/dm6 -U SRR1187947/Trimmed_Reads/SRR1187947_trimmed.fq.gz -S SRR1187947_mapped_verysensitive.sam
now=$(date) 
echo $date > sensitive_end 

now=$(date) 
echo $date > sensitive_local_start 
bowtie2 --very-sensitive-local -a -p 4 -x /Volumes/BlytheLab_Files/HTSeq/Bowtie_Indices/dm6 -U SRR1187947/Trimmed_Reads/SRR1187947_trimmed.fq.gz -S SRR1187947_mapped_verysensitive_local.sam

now=$(date)
echo $date > sensitive_local_end 

