#!/bin/sh

#This script is an automated method to search for Ping Pong Signatures

#conda activate read_search


REFERENCEGFF=/Users/ryan/Documents/GitHub/sm_RNA_seq/sat_dna/wei2020/dmel_scaffold2_plus0310_rm_for_htseq.gff3
bam=/Users/ryan/Documents/GitHub/sm_RNA_seq/sat_dna/mappings/SRR1187947_mapped_verysensitive_local_wei2020_bestreads.mapped.sorted.bam


while IFS= read -r "feature" || [ -n "$feature" ] ; do 

echo $feature

mkdir $feature 

cd $feature 

#Get a bed for the feature 

cat $REFERENCEGFF | grep -F "$feature" | cut -f1 -f4 -f5 -f7 > $feature-region.bed

#Subset the Sam File for the feature 

samtools view -h -L $feature-region.bed $bam > $feature.sam

#Run Signature Analysis

python3 /Users/ryan/Documents/GitHub/sm_RNA_seq/sat_dna/ping_pong/Signature.py $feature.sam 16 32 1 30 $feature.piProfile.txt


#Move the results file 
mv $feature.piProfile.txt ../rsp_piProfiles

#Clean up the extra data
#rm *.sam 
#rm *.bed

#Head back to the starting dir
cd .. ;

#rm -r $feature

done < /Users/ryan/Documents/GitHub/sm_RNA_seq/sat_dna/ping_pong/bp_SAT_elements.txt




#bam=/Users/ryan/Documents/GitHub/sm_RNA_seq/sat_dna/SRR1187947_mapped_verysensitive_local_wei2020_bestreads.mapped.sorted.bam

#REFERENCEGFF=/Users/ryan/Documents/GitHub/sm_RNA_seq/sat_dna/wei2020/dmel_scaffold2_plus0310_rm_for_htseq.gff3



#Get the Regions of interest from the GFF3

#cat $REFERENCEGFF | egrep 'GAGACA' | cut -f1 -f4 -f5 -f7 > $featureregion.bed


#REFERENCEGFF=/Users/ryan/Documents/GitHub/sm_RNA_seq/sat_dna/wei2020/dmel_scaffold2_plus0310_rm_for_htseq.gff3


#samtools view -h -L $bed $bam > $feature.sam



#python3 Signature.py $feature.sam 16 32 1 30 $feature.piProfile.txt

#mv
