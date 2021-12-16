Mirdeep is a software program that identifies novel and known miRNA sequences present in sequencing data. 

This contains code I used to try and run the program. The results were a bit confusing since one a few miRNAs were identified. But I didn't end up trouble shooting further. 

To run mirDeep 

1) Activate the mirDeep conda enviornment 
    ````
    conda activate mirdeep 

    ````
2) Convert your fastq file to a fasta file (gets rid of quality scores) using `fastq_to_fasta.sh` or the following code
    ````

    #Unzip the file if needed 
    gunzip -c $fastq > file

    #Remove the quality scores
    awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' my_file.fq > my_file.fasta

    #Remove the leading whitespace

    perl -i -pe 's|[ \t]||g' my_file.fasta

    ````
3) Convert the BAM file to an ARF file (a special mirDeep file) us the code in `make_arf.sh`
    ```` 
    conda activate mirdeep
    
    #Convert BAM to SAM with the header
    samtools view -h my_file.bam > my_file.sam

    bwa_sam_converter.pl -i my_file.sam -c -o reads_collapsed.fa -a reads_collapsed_vs_genome.arf

    ````
    This will produce a fasta file and arf file that mirdeep can take as an input
4) Run mirdeep

    There are several different variations on how you can actually run mirDeep (with different inputs and different parts of the program). I ran it with the wrapper function `miRDeep2.pl` using this command: 
    ````
    reads=/Users/ryan/Documents/GitHub/sm_RNA_seq/miRNA/reads_collapsed.fa #Removed whitespace
    reference=/Volumes/BlytheLab_Files/HTSeq/Ryan/20210920_smRNAseq_pilot/SRR1187947/Reference/dm6.fa #Removed whitespace
    arf=/Users/ryan/Documents/GitHub/sm_RNA_seq/miRNA/reads_collapsed_vs_genome.arf
    miRNAS=/Users/ryan/Documents/GitHub/sm_RNA_seq/miRNA/mirbase_dme_v22_1.fa

    conda activate mirdeep

    miRDeep2.pl $reads $reference $arf $miRNAS none none -t Fly 
    ````


