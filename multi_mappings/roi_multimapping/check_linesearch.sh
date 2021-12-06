
#Script to if the filterline process ran correctly
for linesfile in /Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/42AB_test/binned_lines2/*; do 
cd $linesfile 

name=$(echo $linesfile | rev | cut -d"/" -f 1 | rev | cut -d"_" -f 1 )

echo $name

cd beds 

lines=$(ls | wc -l)

echo $lines 


done > /Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/42AB_test/linesearchlog.txt


#cat $name $lines /Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/42AB_test/binned_beds/line_searchcheck.txt 