#!/bin/sh


# mkdir GAGA_Repeats 
# cd GAGA_Repeats 
# mkdir raw_data
# cd .. 
# cd SRR1187947/raw_data

# cp SRR1187947.fastq.gz /Users/ryan/Documents/GitHub/sm_RNA_seq/GAGA_Repeats/raw_data

#cd /Users/ryan/Documents/GitHub/sm_RNA_seq/GAGA_Repeats/raw_data




#Unzip the file
#SRR1187947.fastq.gz | head -n 1000 > SRR1187947_test.fastq 




#Get SRA data

out_dir=/Users/ryan/Documents/GitHub/sm_RNA_seq/GAGA_Repeats/raw_data

fasterq-dump -e 4 -O $out_dir -p SRR1187947

#Zip it up 

#gzip $out_dir/*.fastq

#Count the fastq lines
echo $(cat SRR1187947.fastq|wc -l)/4|bc

grep -A 2 -B 1 'AAGAGAAGAGAAGAG' SRR1187947.fastq | sed '/^--$/d' > SRR1187947_AAGAGAAGAGAAGAG.fastq

echo $(cat SRR1187947_AAGAGAAGAGAAGAG.fastq|wc -l)/4|bc
