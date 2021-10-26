#!/bin/sh

fasta=/Users/ryan/Documents/GitHub/sm_RNA_seq/sat_dna/chang2019/dmel_scaffold2_V5.fasta


conda activate mappers


bowtie2-build $fasta chang2019heterochrom

cp chang2019heterochrom*.bt2 /Volumes/BlytheLab_Files/HTSeq/Bowtie_Indices