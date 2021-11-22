for linesfile in /Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/histone_cluster/binned_lines/*; do 

name=$(echo $linesfile | rev | cut -d"/" -f 1 | rev | cut -d"_" -f 1 )

cd $linesfile 

cd beds 

rm *.bed

cd .. # leave the beds file 
rm -r beds

cd .. # leave the lines file 
done 