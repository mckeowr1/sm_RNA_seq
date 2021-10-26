#!/bin/sh


wang_gtf=/Users/ryan/Documents/GitHub/sm_RNA_seq/sat_dna/wei2020/dmel_scaffold2_plus0310_rm_for_htseq.gff3

#AAGAG Repeats 
cat $wang_gtf | grep -F '(AAGAG)n' > aagag_repeats.gff3

