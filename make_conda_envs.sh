#!/bin/sh
#Make Conda EnV

#Put Channels In the Right order for Bioconda
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda create -n mappers bowtie2

conda activate mappers      