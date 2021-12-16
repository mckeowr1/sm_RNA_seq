##################
### Map to DM6 ### 
##################

#Map a 50mer of AAGAG
bowtie2 --very-sensitive-local -a -p 4 -x /Volumes/BlytheLab_Files/HTSeq/Bowtie_Indices/dm6 -U /Users/ryan/Documents/GitHub/sm_RNA_seq/GAGA_Repeats/raw_data/gaga50.fastq 2> GAGA_Repeats.log.txt |
samtools view - > gaga_verysensitive_local.mapped.sam

#Count number of mappings
wc -l gaga_verysensitive_local.mapped.sam
#1358

#Map a 50mer of AACAC 
bowtie2 --very-sensitive-local -a -p 4 -x /Volumes/BlytheLab_Files/HTSeq/Bowtie_Indices/dm6 -U /Users/ryan/Documents/GitHub/sm_RNA_seq/GAGA_Repeats/raw_data/aacac50.fastq 2> AACAC_Repeats.log.txt |
samtools view - > aacac_verysensitive_local.mapped.sam


wc -l aacac_verysensitive_local.mapped.sam
#59

#Map a 50mer of AATAT
bowtie2 --very-sensitive-local -a -p 4 -x /Volumes/BlytheLab_Files/HTSeq/Bowtie_Indices/dm6 -U /Users/ryan/Documents/GitHub/sm_RNA_seq/GAGA_Repeats/raw_data/aatat50.fastq 2> aatat_Repeats.log.txt |
samtools view - > aatat_verysensitive_local.mapped.sam

#2664

#########################
### Map to Larracuente ### 
#########################

index=/Volumes/BlytheLab_Files/HTSeq/Bowtie_Indices/chang2019heterochrom

#Map a 50mer of AAGAG
bowtie2 --very-sensitive-local -a -p 4 -x $index -U /Users/ryan/Documents/GitHub/sm_RNA_seq/GAGA_Repeats/raw_data/gaga50.fastq 2> GAGA_Repeats_heterochrom.log.txt |
samtools view -bS > gaga_verysensitive_local_heterochrom.mapped.sam

#Map a 50mer of AACAC
bowtie2 --very-sensitive-local -a -p 4 -x $index -U /Users/ryan/Documents/GitHub/sm_RNA_seq/GAGA_Repeats/raw_data/aacac50.fastq 2> aacac_Repeats_heterochrom.log.txt |
samtools view -bS > aacac_verysensitive_local_heterochrom.mapped.sam

#Map a 50mer of AATAT
bowtie2 --very-sensitive-local -a -p 4 -x $index -U /Users/ryan/Documents/GitHub/sm_RNA_seq/GAGA_Repeats/raw_data/aatat50.fastq 2> aatat_Repeats_heterochrom.log.txt |
samtools view -bS > aatat_verysensitive_local_heterochrom.mapped.sam
