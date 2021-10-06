#!/bin/sh

#Define the search beds 

bedaa=/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/search_beds/SRR1187947_mapped_verysensitive_local_sortedaa
bedab=/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/search_beds/SRR1187947_mapped_verysensitive_local_sortedab
bedac=/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/search_beds/SRR1187947_mapped_verysensitive_local_sortedac
bedad=/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/search_beds/SRR1187947_mapped_verysensitive_local_sortedad
bedae=/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/search_beds/SRR1187947_mapped_verysensitive_local_sortedae
bedaf=/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/search_beds/SRR1187947_mapped_verysensitive_local_sortedaf
bedag=/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/search_beds/SRR1187947_mapped_verysensitive_local_sortedag
bedah=/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/search_beds/SRR1187947_mapped_verysensitive_local_sortedah

cd /Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/42AB_test/binned_lines

for linesfile in /Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/42AB_test/binned_lines/*; do 
cd $linesfile 
aalines=$linesfile/*_aa.txt
ablines=$linesfile/*_ab.txt 
aclines=$linesfile/*_ac.txt
adlines=$linesfile/*_ad.txt 
aelines=$linesfile/*_ae.txt 
aflines=$linesfile/*_af.txt 
aglines=$linesfile/*_ag.txt 
ahlines=$linesfile/*_ah.txt

#name=$(echo $linesfile | cut -d_ -f1)

mkdir beds

/Users/ryan/filterline/filterline $aalines $bedaa > beds/aa.bed & 
/Users/ryan/filterline/filterline $ablines $bedab > beds/ab.bed & 
/Users/ryan/filterline/filterline $aclines $bedac > beds/ac.bed & 
/Users/ryan/filterline/filterline $adlines $bedad > beds/ad.bed & 
/Users/ryan/filterline/filterline $aelines $bedae > beds/ae.bed & 
/Users/ryan/filterline/filterline $aflines $bedaf > beds/af.bed & 
/Users/ryan/filterline/filterline $aglines $bedag > beds/ag.bed & 
/Users/ryan/filterline/filterline $ahlines $bedah > beds/ah.bed

cat a*.bed > $name.bed 

mv $name.bed 42AB_beds


cd .. 
done




# bedaa=/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/42AB/SRR1187947_mapped_verysensitive_local_sortedaa
# bedac=/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/42AB/SRR1187947_mapped_verysensitive_local_sortedac 

# /Users/ryan/filterline/filterline /Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/42AB/bin10_lines_aa.txt $bedaa > bin10_bed_aa & 
# /Users/ryan/filterline/filterline /Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/42AB/bin10_lines_ac.txt $bedac > bin10_bed_ac

