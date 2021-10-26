#!/bin/sh

#Define the search beds 

bedaa=/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/search_beds/SRR1187947_mapped_verysensitive_local_sortedaa.bed
bedab=/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/search_beds/SRR1187947_mapped_verysensitive_local_sortedab.bed
bedac=/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/search_beds/SRR1187947_mapped_verysensitive_local_sortedac.bed
bedad=/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/search_beds/SRR1187947_mapped_verysensitive_local_sortedad.bed
bedae=/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/search_beds/SRR1187947_mapped_verysensitive_local_sortedae.bed
bedaf=/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/search_beds/SRR1187947_mapped_verysensitive_local_sortedaf.bed
bedag=/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/search_beds/SRR1187947_mapped_verysensitive_local_sortedag.bed
bedah=/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/search_beds/SRR1187947_mapped_verysensitive_local_sortedah.bed

cd /Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/histone_cluster/binned_lines

for linesfile in /Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/histone_cluster/binned_lines4/*; do 
cd $linesfile 
aalines=$linesfile/*_aa.txt
ablines=$linesfile/*_ab.txt 
aclines=$linesfile/*_ac.txt
adlines=$linesfile/*_ad.txt 
aelines=$linesfile/*_ae.txt 
aflines=$linesfile/*_af.txt 
aglines=$linesfile/*_ag.txt 
ahlines=$linesfile/*_ah.txt

name=$(echo $linesfile | rev | cut -d"/" -f 1 | rev | cut -d"_" -f 1 )

echo $name
mkdir beds

/Users/ryan/filterline/filterline $aalines $bedaa > beds/aa.bed & 
/Users/ryan/filterline/filterline $ablines $bedab > beds/ab.bed & 
/Users/ryan/filterline/filterline $aclines $bedac > beds/ac.bed & 
/Users/ryan/filterline/filterline $adlines $bedad > beds/ad.bed & 
/Users/ryan/filterline/filterline $aelines $bedae > beds/ae.bed & 
/Users/ryan/filterline/filterline $aflines $bedaf > beds/af.bed & 
/Users/ryan/filterline/filterline $aglines $bedag > beds/ag.bed & 
/Users/ryan/filterline/filterline $ahlines $bedah > beds/ah.bed

#cd beds

#cat a*.bed > $name.bed 

#cp $name.bed /Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/42AB_test/42AB_beds

#cd .. #Get out of bed directory
cd .. #Get out of Lines directory

done






