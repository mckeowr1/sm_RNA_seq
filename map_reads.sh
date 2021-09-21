#!/bin/sh

## Parse Command Line Args ## 

basedir=$1 #Project Directory 
index=$2 #Path to Index Files in NSA
filename=$3 #Trimmed reads file

cd $basedir


#Make a Dir for the Mapping 
mkdir Mapped_Reads

#Define sub-dir structure
trimdir=$base_dir/Trimmed_Reads
mapdir=$base_dir/Mapped_Reads


#Make a variable to name files with  
shortname=$(echo $filename |cut -d_ -f1-2)

#This is where the problem is!
bowtie2 -p 4 -X 2000 -x $index -1 $basedir/$trimdir/$filename -2 $basedir/$trimdir/$filename2 2> $basedir/$mapdir/$shortname.bowtie.mapping.log.txt |  \ #Create a log of reads that map and reads that don't
samtools view -bS - > $basedir/$mapdir/$shortname.mapped.bam #Create a bam file
samtools sort -@ 4 -m 2G -T $basedir/$mapdir/temp -o $basedir/$mapdir/temp1.bam $basedir/$mapdir/$shortname.mapped.bam



# #This is where things have to change!

# java -jar /Applications/picard.jar CleanSam \ #Clean the Sam?
# I=$basedir/$mapdir/temp1.bam \
# O=$basedir/$mapdir/temp2.bam
# java -jar /Applications/picard.jar MarkDuplicates \
# I=$basedir/$mapdir/temp2.bam \
# O=$basedir/$mapdir/$shortname.mapped.md.bam \
# M=$basedir/$mapdir/$shortname.dup.metrics.txt
# rm $basedir/$mapdir/temp1.bam $basedir/$mapdir/temp2.bam $basedir/$mapdir/$shortname.mapped.bam
# samtools index $basedir/$mapdir/$shortname.mapped.md.bam