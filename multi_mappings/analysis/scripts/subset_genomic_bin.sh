#!/bin/zsh
bed_dir=/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/histone_cluster/histone_cluster_beds



while IFS= read -r "binned_bed" || [ -n "$binned_bed" ] ; do 


cd bed_dir

bed=$(ls | grep $binned_bed)

# zshell way to get path for the file
mypath=${histone_bed:a}





bedtools intersect -a bin260.bed -b test_region.bed   


echo $binned_bed 

done < test_beds.txt



#bedtools intersect -a bin260.bed -b test_region.bed      