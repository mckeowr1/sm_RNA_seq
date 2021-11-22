#!/bin/zsh

## This script is used to subset the reads from a particular ROI that also map to another ROI
## Ex) Reads that map to Bin10 of the histone cluster and also map to a particular Gene ex)Luna 
## 1st all reads that map to the particular genomic loci of interest are pulled
## Then the names of the reads for that loci that also map to the ROI are pulled by subsetting the ROI_bin_bed using bedtools intersect
## Next picard tools is used to filter the genomic loci bed to just reads that map to both locations 
## We can then vizualize the reads using IGV (this requires sorting and indexing with samtools)
## Or we can look for regulatory signatures like Ping Pong


#Conda activate mappers

# Get All reads that map to genomic loci
bedtools intersect -a /Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/files/SRR1187947_mapped_verysensitive_local.mapped.bam -b FBgn0040765_roi.bed > FBgn0040765.sam 



#Make a list of reads that map to ROI and genomic loci
bedtools intersect -a bin263.bed -b fbgn_0027339_roi.bed  > fbgn_0027339_histone_reads.bed #Get reads that map to the ROI Bin and are assoicated with another genomic loci ex) Luna


#Process the reads into a READ_LIST_FILE 

Rscript region_bams.R fbgn_0027339_histone_reads.bed


java -jar /Applications/picard.jar FilterSamReads I=FBgn0040765.bam O=FBgn0040765.bam_histone.bam READ_LIST_FILE=FBgn0040765_histone_reads.txt FILTER=includeReadList


python3 /Users/ryan/Documents/GitHub/sm_RNA_seq/sat_dna/ping_pong/Signature.py FBgn0040765.bam 16 32 1 30 fbn_0040765_histone.piProfile.txt



#bedtools intersect -a bin260.bed -b test_region.bed      