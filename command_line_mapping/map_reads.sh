#!/bin/sh

## Parse Command Line Args ## 

basedir=$1 #Project Directory 
index=$2 #Path to Index Files in NSA
filename=$3 #Trimmed reads FASTQ file


cd $basedir


#Make a Dir for the Mapping 
mkdir Mapped_Reads

#Define sub-dir structure
trimdir=Trimmed_Reads
mapdir=Mapped_Reads


#Make a variable to name files with  
shortname=$(echo $filename |cut -d_ -f1-2)
echo $shortname
#conda activate mappers 

##Bowtie2 Mapping 
#-X maximum insert size for paired end
#-x bowtie index
#-q query file 

bowtie2 -p 4 -x $index -U $trimdir/$filename 2> $mapdir/$shortname.bowtie.mapping.log.txt |
samtools view -bS - > $mapdir/$shortname.mapped.bam
samtools sort -@ 4 -m 2G -T $mapdir/temp -o $mapdir/temp1.bam $mapdir/$shortname.mapped.bam


##Let's Clean things up a bit ! 
#Clean Sam Fixes the header
java -jar /Applications/picard.jar CleanSam \
I=$mapdir/temp1.bam \
O=$mapdir/temp2.bam

#MarkDuplicates
java -jar /Applications/picard.jar MarkDuplicates \
I=$mapdir/temp2.bam \
O=$mapdir/$shortname.mapped.md.bam \
M=$mapdir/$shortname.dup.metrics.txt


#Get Rid of the Temporary Files
rm $mapdir/temp1.bam $mapdir/temp2.bam $mapdir/$shortname.mapped.bam

samtools index $mapdir/$shortname.mapped.md.bam


#conda deactivate





