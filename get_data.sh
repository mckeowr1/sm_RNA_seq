#!/bin/sh

#Get SRA data

out_dir=/Volumes/BlytheLab_Files/HTSeq/Ryan/20210920_smRNAseq_pilot/SRR1187947/raw_data

fasterq-dump -e 4 -O $out_dir -p SRR1187947

#Zip it up 

gzip $out_dir/*.fastq