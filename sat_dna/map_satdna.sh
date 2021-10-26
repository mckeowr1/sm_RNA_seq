#!/bin/sh

conda activate mappers

trimfastq=/Users/ryan/Documents/GitHub/sm_RNA_seq/SRR1187947/Trimmed_Reads/SRR1187947_trimmed.fq
index=/Volumes/BlytheLab_Files/HTSeq/Bowtie_Indices/chang2019heterochrom

#Get all reads
bowtie2 --very-sensitive-local -a -p 4 -x $index -U $trimfastq 2> SRR1187947_mapped_verysensitive_local_repBase_allreads.log.txt |
samtools view -bS - > SRR1187947_mapped_verysensitive_local_repBase_allreads.mapped.bam

#Random selection of best mapping read 
bowtie2 --very-sensitive-local -p 4 -x $index -U $trimfastq 2> SRR1187947_mapped_verysensitive_local_repBase_bestreads.log.txt |
samtools view -bS - > SRR1187947_mapped_verysensitive_local_repBase_bestreads.mapped.bam 

samtools sort -@ 8 -m 2G -T temp -o temp1.bam SRR1187947_mapped_verysensitive_local_wei2020_bestreads.mapped.bam


java -jar /Applications/picard.jar CleanSam \
I=temp1.bam \
O=temp2.bam


#Skipping mark duplicates 

#Index the bam


