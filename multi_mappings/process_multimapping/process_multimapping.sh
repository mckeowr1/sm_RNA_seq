#!/bin/bash
#SBATCH -J make        ## Name of job
#SBATCH -A b1042               ## Allocation
#SBATCH -p genomics               ## Queue
#SBATCH -t 5:00:00             ## Walltime/duration of the job
#SBATCH --mem=32G           ## Total memory in GB needed for a job.
#SBATCH --cpus-per-task=4           ## Number of processors for your task

bam=path/to/your/multimapped_bam

name=$(echo $bam | rev | cut -d"/" -f 1 | rev )

module load bedtools

bedtoosl bamtobed -tag AS -i $bam > $name.bed

bedtools sort -i $name.bed

