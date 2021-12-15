#!/bin/bash
#SBATCH -J clean_sam       ## Name of job
#SBATCH -A b1042               ## Allocation
#SBATCH -p genomics               ## Queue
#SBATCH -t 5:00:00             ## Walltime/duration of the job
#SBATCH --mem=30G           ## Total memory in GB needed for a job.
#SBATCH --cpus-per-task=4           ## Number of processors for your task

bed=/path/to/your/multimapped_bed

split --lines=100000000 $bed

