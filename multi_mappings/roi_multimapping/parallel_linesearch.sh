#!/bin/zsh

#Define some paths
search_path=/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/search_beds
lines_dir=/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/histone_cluster/binned_lines #This should be the set directory you want to run

#Define the search beds 

bedaa="$search_path"/*aa.bed
bedab="$search_path"/*ab.bed
bedac="$search_path"/*ac.bed
bedad="$search_path"/*ad.bed
bedae="$search_path"/*ae.bed
bedaf="$search_path"/*af.bed
bedag="$search_path"/*ag.bed
bedah="$search_path"/*ah.bed
#Make sure to add a variable if you have more than 8 search beds


cd $lines_dir


for linesfile in "$lines_dir"/*; do 
cd $linesfile 
aalines=$linesfile/*_aa.txt
ablines=$linesfile/*_ab.txt  
aclines=$linesfile/*_ac.txt
adlines=$linesfile/*_ad.txt 
aelines=$linesfile/*_ae.txt 
aflines=$linesfile/*_af.txt 
aglines=$linesfile/*_ag.txt 
ahlines=$linesfile/*_ah.txt
#Make sure to add a variable if you have more than 8 search beds


name=$(echo $linesfile | rev | cut -d"/" -f 1 | rev | cut -d"_" -f 1 )

echo $name
mkdir beds


#Make sure to edit the path to your filterline program
/Users/ryan/filterline/filterline $aalines $bedaa > beds/aa.bed & 
/Users/ryan/filterline/filterline $ablines $bedab > beds/ab.bed & 
/Users/ryan/filterline/filterline $aclines $bedac > beds/ac.bed & 
/Users/ryan/filterline/filterline $adlines $bedad > beds/ad.bed & 
/Users/ryan/filterline/filterline $aelines $bedae > beds/ae.bed & 
/Users/ryan/filterline/filterline $aflines $bedaf > beds/af.bed & 
/Users/ryan/filterline/filterline $aglines $bedag > beds/ag.bed & 
/Users/ryan/filterline/filterline $ahlines $bedah > beds/ah.bed
#Make sure to add another search program if you have more than 8 search beds


#cd beds

#cat a*.bed > $name.bed 

#cp $name.bed /Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/42AB_test/42AB_beds

#cd .. #Get out of bed directory
cd .. #Get out of Lines directory

done






