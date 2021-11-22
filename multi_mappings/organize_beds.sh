for linesfile in /Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/histone_cluster/binned_lines/*; do 

name=$(echo $linesfile | rev | cut -d"/" -f 1 | rev | cut -d"_" -f 1 )

cd $linesfile 

cd beds 

cat a*.bed > $name.bed 

cp $name.bed /Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/histone_cluster/histone_cluster_beds

cd .. 
cd ..
done 