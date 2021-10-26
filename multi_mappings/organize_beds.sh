for linesfile in /Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/42AB_test/binned_lines4/*; do 

name=$(echo $linesfile | rev | cut -d"/" -f 1 | rev | cut -d"_" -f 1 )

cd $linesfile 

cd beds 

cat a*.bed > $name.bed 

cp $name.bed /Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/42AB_test/42AB_beds

cd .. 
cd ..
done 