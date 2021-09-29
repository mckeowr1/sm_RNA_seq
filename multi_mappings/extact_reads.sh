#!/bin/sh

for filename in /Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/42AB_Lines/*.txt ;
do
echo $filename

/Users/ryan/filterline/filterline $filename files/SRR1187947_mapped_verysensitive_local.mapped.bed > $filename.bed


done

#/Users/ryan/filterline/filterline 42AB_Lines/SRR1187947.99_lines.txt files/SRR1187947_mapped_verysensitive_local.mapped.bed

