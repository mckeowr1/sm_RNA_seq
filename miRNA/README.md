Mirdeep is a software program that identifies novel and known miRNA sequences present in sequencing data. 

This contains code I used to try and run the program. The results were a bit confusing since one a few miRNAs were identified. But I didn't end up trouble shooting further. 

To run mirDeep 

1) Activate the mirDeep conda enviornment 
    ````
    conda activate mirdeep 

    ````
2) Convert your fastq file to a fasta file (gets rid of quality scores) using `fastq_to_fasta.sh` 
````
#Unzip the file if needed 
gunzip -c $fastq > file

#Remove the quality scores
awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' my_file.fq > my_file.fasta

#Remove the leading whitespace

perl -i -pe 's|[ \t]||g' my_file.fasta


````
