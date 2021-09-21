# sm_RNA_seq
Blythe Lab GAGA smRNA seq 

# Data
Raw-Data: /Volumes/BlytheLab_Files/HTSeq/Ryan/20210920_smRNAseq_pilot/SRR1187947/raw_data
Notes: 
single read 


# Workflow 
1) Get Data from SRA 
2) Trim the reads
3) Map the reads to Genome of Interest



## Trim Reads Dir Structure

project
    Raw_Data
    Trimmed_Reads

# Running Trim Reads Executable 

`bash trim_reads.sh $projectdir`

# Dir Structure Should Look Like this 
project
    Raw_Data
    Trimmed_Reads
        SRR1187947_test.fastq.gz_trimming_report.txt
        SRR1187947_test_trimmed.fq.gz
        SRR1187947_test_trimmed_fastqc.html #HTML Report with Trimming statst
        SRR1187947_test_trimmed_fastqc.zip


# Do we want any checks from the Trimming


## Mapping Reads 

Conda Env
`bowtie2 2.4.4`

`conda activate mappers`
`bash map_reads.sh project_dir path/to/bt_index FASTQ_file` 


project
    Raw_Data
    Trimmed_Reads
    Mapped_Reads

`bash map_reads.sh $project_dir $ref_index` #Will need some way to handle non DM6 index

path to DM6 Index:  `/Volumes/BlytheLab_Files/HTSeq/Bowtie_Indices/dm6` 



#Don't mark duplicates w/ picard tools 


# Need to build Bowtie Indices for Larracuente 2017


bowtie build manual 




path to RepBase Genome `afp://blythelabnas.mbs.northwestern.edu/BlytheLab_Files/HTSeq/Genomes/Drosophila_Repbase_Khost_Larracuente_2017` 






