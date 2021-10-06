#!/bin/sh 

bedtools bamtobed 


bedtools bamtobed -tag AS -i SRR1187947_mapped_verysensitive_local.mapped.bam > SRR1187947_mapped_verysensitive_local.mapped.bed


bedtools sort -i test.bed


#Get a Region of interest (Cluster42AB) from the bed 

bedtools intersect -a test.bed -b region.bed