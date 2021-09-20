#!/bin/sh

basedir=$1

cd basedir

mkdir Trimmed_Reads 


rawdir=Raw_Data 
trimdir=Trimmed_Reads


#Activate trimgalor env
source /Applications/TrimGaloreEnv/bin/activate


#Check for paired end 
#gzcat SRR1187947.fastq.gz | wc -l 





### Trim the Reads ### 
#######################


trim_galore -o $basedir/$trimdir --gzip --fastqc $file 


#for file in $basedir/$rawdir/*.fastq.gz
#do 
#filename=$(basename "$file")
#echo "path to the first read"
#echo $basedir/$rawdir/$filename
#echo "base name for the first read"
#echo $filename
#filename2=$(echo $filename |cut -d_ -f1-3|awk '{print $1"_R2_001.fastq.gz"}')
#echo "path to the second read"
#echo $basedir/$rawdir/$filename2
#echo "base name for the second read"
#echo $filename2
#echo
#echo
#trim_galore -o $basedir/$trimdir --gzip --fastqc --paired $basedir/$rawdir/$filename $basedir/$rawdir/$filename2
#done


deactivate #Get Rid of Trim Galore ENV








