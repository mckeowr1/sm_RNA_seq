#!/bin/zsh

bedtools intersect -a bin1139.bed -b h4_roi.bed  > h4_histone_reads.bed 

bedtools intersect -a /Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/files/SRR1187947_mapped_verysensitive_local.mapped.bam -b h4_roi.bed > h4.bam



java -jar /Applications/picard.jar FilterSamReads I=h4.bam O=h4_histone.bam READ_LIST_FILE=h4_histone_reads.txt FILTER=includeReadList


python3 /Users/ryan/Documents/GitHub/sm_RNA_seq/sat_dna/ping_pong/Signature.py fbgn_0027339.sam 16 32 1 30 fbn_0027339_histone.piProfile.txt



#bedtools intersect -a bin260.bed -b test_region.bed      



## Code to loop through for multiple files

#bed_dir=/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/histone_cluster/histone_cluster_beds
#while IFS= read -r "binned_bed" || [ -n "$binned_bed" ] ; do 
#cd bed_dir
#bed=$(ls | grep $binned_bed)
# zshell way to get path for the file
#mypath=${histone_bed:a}
#bedtools intersect -a bin260.bed -b test_region.bed   
#echo $binned_bed 

#done < test_beds.txt



#bedtools intersect -a bin260.bed -b test_region.bed      