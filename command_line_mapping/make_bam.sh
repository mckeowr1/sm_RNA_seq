#!/bin/sh 


#This takes the output of bowtie and makes a bam out of it to read into R 
#This is the same as map_reads.sh but is bowtie has already been performed


basedir=$1 #Path to Project Dir
sam=$2 #File that will be BAMed

shortname=$(echo $sam| cut -d. -f1 )

cd $basedir

mkdir BAMS

bamdir=BAMS


samtools view -bS - $sam -> $bamdir/$shortname.mapped.bam
samtools sort  -@ 4 -m 2G -T $bamdir/temp -o $bamdir/temp1.bam $bamdir/$shortname.mapped.bam


##Let's Clean things up a bit ! 
#Clean Sam Fixes the header
java -jar /Applications/picard.jar CleanSam \
I=$bamdir/temp1.bam \
O=$bamdir/temp2.bam


java -jar /Applications/picard.jar MarkDuplicates \
I=$bamdir/temp2.bam \
O=$bamdir/$shortname.mapped.md.bam \
M=$bamdir/$shortname.dup.metrics.txt

rm $bamdir/temp1.bam $bamdir/temp2.bam $bamdir/$shortname.mapped.bam

samtools index $bamdir/$shortname.mapped.md.bam

