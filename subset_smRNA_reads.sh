#!/bin/sh 

#Take an input bam & bed file and subsets the reads that map to a Regiong of interest 

basedir=$1 #A project directory that contains a Mapped Reads file with a BAM in it (Not marked for duplicates)
bed=$2 #A bed file that contains the regions of interest

#Get the name of the bead file
roiname=$(echo $bed | cut -d. -f1)

mkdir $roiname #Make a directory for the region of interest 

cd $basedir

#If we are looking for a particular strand this may not work! 
samtools view -h -L $bed SRR1187947_mapped_verysensitive_local.mapped.bam > $roiname.sam

#Trying with bedtools 
#bedtools intersect -abam SRR1187947_mapped_verysensitive_local.mapped.bam -b 42AB.bed -ubam

#Check if the file is empty?
head_size=$(samtools view -H SRR1187947_mapped_verysensitive_local.mapped.bam | wc -l)
echo "There are" $head_size "lines in the header"

roi_bam_size=$(samtools view -H SRR1187947_mapped_verysensitive_local.mapped.bam | wc -l)
##Where else do reads map## 

#Subset the reads that have equally as good mappings
    #Subset bu flags or by ASi = XS:I

#Filter By Flag 256 Not Primary Mapping 
#To get the primary mapping can pull from another sam file where are the reads aren't listed

#Potentially implement a quality matching 

#Get all Non-Primary Mapping Reads in the region
samtools view -h -f 256 $roiname.sam > $roiname.multi_mapping.sam

#Spit out a list of read names
samtools view $roiname.multi_mapping.sam| cut -f1 | sort | uniq -u > $roiname-multimapping_read_names.txt



## Loop through Reads to Pull them from Original BAM ## 

# for read_name
# do 
# grep '$read_name' bam > $read_name.bam #A bam file the other places that read maps

# grep -w 'SRR1187947.7788'  

#This step may need to get alot faster! - 
samtools view SRR1187947_mapped_verysensitive_local.mapped.bam | grep -w SRR1187947.57 > SRR1187947_read.sam

#Picard 

# java -jar /Applications/picard.jar MarkDuplicates \
# I=$mapdir/temp2.bam \
# O=$mapdir/$shortname.mapped.md.bam \
# M=$mapdir/$shortname.dup.metrics.txt

java -jar /Applications/picard.jar FilterSamReads \
I=SRR1187947_mapped_verysensitive_local.mapped.bam \
O=output.bam \
READ_LIST_FILE=read.txt \
FILTER=includeReadList

##Take the read BAM and Intersect with a GFF ## 

#Output a GFF file with regions that the reads in the BAM overlap with 


bedtools intersect -a rev2.gff -b SRR1187947.57.bam -wa 

head -n 28100000