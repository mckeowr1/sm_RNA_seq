# Input Dir Structure
/42AB_Test
├── bed_roi
│   ├── 42AB.bed
├── search_beds
├── binned_lines
├── binned_beds
└── binned_readnames



# Output Dir Structure


/42AB_Test
├── bed_roi
│   ├── 42AB.bed
├── search_beds
│   ├── SRR1187947_mapped_verysensitive_local.mapped_sorted.bedindex.tsv
    ├── SRR1187947_mapped_verysensitive_local.mapped_sorted.bed 
    ├── SRR1187947_mapped_verysensitive_local.mapped_sortedaa.bed
    ├── SRR1187947_mapped_verysensitive_local.mapped_sortedab.bed
└── binned_readnames
│   ├── bin1_readnames.txt
│   ├── bin2_readnames.txt
│   ├── ...    
├── binned_lines
    ├── bin1_lines
        ├──beds 
            ├──aa.bed
            ├──ab.bed
            ├──ac.bed
            ├──ad.bed
            ├──ae.bed
            ├──af.bed
            ├──ag.bed
            ├──ah.bed
        ├──bin1_aa.txt
        ├──bin1_ab.txt
        ├──bin1_ac.txt
        ├──bin1_ad.txt
        ├──bin1_ae.txt
        ├──bin1_af.txt
        ├──bin1_ag.txt
        ├──bin1_ah.txt
├── design
    └── 20191119_design.csv


# Scripts

### multi_map_matrix.R
Bins the region of interest and produces list of readnames in each bin

Inputs: 
A bed region file that has been created from the mapped bams
Chrom, Start, Stop of Region

Outputs: 
A directory with readnames text files 

### find_reads.R 
Pulls the line location of a list of reads 

Inputs: 
Sorted Bed Index 
Directory with readnames text files

Output:
Lines Directory
