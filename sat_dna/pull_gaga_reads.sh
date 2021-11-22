#!/bin/sh



#Pull GAGA orginating reads

#conda activate read_search 


#roi_bed=/Users/ryan/Documents/GitHub/sm_RNA_seq/sat_dna/aagag_repeats.bed
bam=/Users/ryan/Documents/GitHub/sm_RNA_seq/sat_dna/SRR1187947_mapped_verysensitive_local_wei2020_bestreads.mapped.sorted.bam


#bedtools intersect -abam $bam -b $roi_bed 


extract_out=$workingdir/extract_reads/${samplename}

REFERENCEGFF=/Users/ryan/Documents/GitHub/sm_RNA_seq/sat_dna/wei2020/dmel_scaffold2_plus0310_rm_for_htseq.gff3


python3 extract_sequence_by_feature_gff.py -t "GGCGAA" "$bam" "$REFERENCEGFF" "GGCGAA_trimmed.fq"



cat bp_SAT_trimmed.fq | grep -A1 "^[>@]" | sed 's/@/>/g'| grep -v "^--" > "bp_SAT_trimmed.fasta"


#Convert the reads to a SAM File 

conda activate mappers

index=/Volumes/BlytheLab_Files/HTSeq/Bowtie_Indices/chang2019heterochrom
bowtie2 --very-sensitive-local -p 4 -x $index -f bp_SAT_trimmed.fasta 2> bp_SAT_mapped_verysensitive_local_repBase_bestreads.log.txt |
samtools view -bS - > bp_SAT_mapped_verysensitive_local_repBase_bestreads.mapped.bam 

samtools sort -@ 8 -m 2G -T temp -o temp1.bam bp_SAT_mapped_verysensitive_local_repBase_bestreads.mapped.bam

java -jar /Applications/picard.jar CleanSam \
I=temp1.bam \
O=bp_SAT_mapped_verysensitive_local_repBase_bestreads.mapped.clean.bam


samtools index bp_SAT_mapped_verysensitive_local_repBase_bestreads.mapped.clean.bam 


samtools view -h bp_SAT_mapped_verysensitive_local_repBase_bestreads.mapped.clean.bam > bp_SAT_mapped_verysensitive_local_repBase_bestreads.mapped.clean.sam

#Check for the Ping Pong Profile

conda activate read_search

python3 Signature.py bp_SAT_mapped_verysensitive_local_repBase_bestreads.mapped.clean.sam 18 40 1 30 bp_SAT.piProfile.txt



