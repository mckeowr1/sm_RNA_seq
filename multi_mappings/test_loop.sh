#!/bin/sh

mkdir sequence_bams

while IFS= read -r seq; do 
echo $seq > read_name.txt
head read_name.txt

java -jar /Applications/picard.jar FilterSamReads \
I=SRR1187947_mapped_verysensitive_local.mapped.bam \
O=sequence_bams/$seq.bam \
READ_LIST_FILE=read_name.txt \
FILTER=includeReadList

done < test_readnames.txt
