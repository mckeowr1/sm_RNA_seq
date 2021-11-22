#!/bin/sh

conda activate read_search

bam=/Users/ryan/Documents/GitHub/sm_RNA_seq/sat_dna/SRR1187947_mapped_verysensitive_local_wei2020_bestreads.mapped.sorted.bam
REFERENCEGFF=/Users/ryan/Documents/GitHub/sm_RNA_seq/sat_dna/wei2020/dmel_scaffold2_plus0310_rm_for_htseq.gff3
GAGAGFF=/Users/ryan/Documents/GitHub/sm_RNA_seq/sat_dna/aagag_repeats.gff3
editGFF=/Users/ryan/Documents/GitHub/sm_RNA_seq/sat_dna/wei2020/dmel_scaffold2_plus0310_rm_for_htseq.edit.gff3

python3 htseq_bam_count_proportional.py $bam $REFERENCEGFF count.out