reads=/Users/ryan/Documents/GitHub/sm_RNA_seq/miRNA/reads_collapsed.fa #Removed whitespace
reference=/Volumes/BlytheLab_Files/HTSeq/Ryan/20210920_smRNAseq_pilot/SRR1187947/Reference/dm6.fa #Removed whitespace
arf=/Users/ryan/Documents/GitHub/sm_RNA_seq/miRNA/reads_collapsed_vs_genome.arf
miRNAS=/Users/ryan/Documents/GitHub/sm_RNA_seq/miRNA/mirbase_dme_v22_1.fa

conda activate mirdeep

miRDeep2.pl $reads $reference $arf $miRNAS none none -t Fly 