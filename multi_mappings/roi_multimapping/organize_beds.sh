#!/bin/zsh

#Define the line & output directory 
lines_dir=/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/histone_cluster/binned_lines #This should be the set directory you want to run
out_dir=/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/histone_cluster/histone_cluster_beds

for linesfile in "$lines_dir"/*; do 

name=$(echo $linesfile | rev | cut -d"/" -f 1 | rev | cut -d"_" -f 1 )

cd $linesfile 

cd beds 

cat a*.bed > $name.bed 

cp $name.bed $out_dir

cd .. 
cd ..
done 