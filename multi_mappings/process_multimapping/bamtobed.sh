#!/bin/sh 

#conda activate mappers

bam=/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/histone_cluster/histonecluster.sam
outname=histone_cluster/histonecluster_reads.bed

bedtools bamtobed -tag AS -i $bam > $outname

#edtools sort -i test.bed


#Get a Region of interest (Cluster42AB) from the bed 

#bedtools intersect -a test.bed -b region.bed