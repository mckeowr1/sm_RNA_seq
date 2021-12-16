---
title: "Trimming and mapping reads"
author: "Shelby Blythe"
date: "9/2/2020"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, tidy = TRUE)
```

This is a worked example for trimming and mapping reads. You may have to adapt this code depending on the type of experiment you are doing. 

In this example, we have paired-end reads from a MiSeq run. The raw data have been downloaded to the NAS. The MiSeq puts each sample in its own directory, and I have manually pulled out the barcode-split reads into a single folder called "Raw_Data."

```
"/Volumes/BlytheLab_Files/HTSeq/Shelby/200824_cg_seChIP_pilot/Raw_Data"
```

Above is the path to the raw data. One step back in the path is the directory for this analysis, we will be making the following additional subdirectories: "Trimmed_Reads" and "Mapped_Reads"

```{bash, eval = FALSE}
mkdir /Volumes/BlytheLab_Files/HTSeq/Shelby/200824_cg_seChIP_pilot/Trimmed_Reads
mkdir /Volumes/BlytheLab_Files/HTSeq/Shelby/200824_cg_seChIP_pilot/Mapped_Reads
```

What we have in the `Raw_Data` directory is the following:

```{bash}
ls /Volumes/BlytheLab_Files/HTSeq/Shelby/200824_cg_seChIP_pilot/Raw_Data
```

What we want to do is to cycle through these files, but because this is paired-end, we need to feed analysis scripts the file for the first and second read. To do this, we employ a 'for-loop' in `bash` using the names of the first reads, and then use a combination of `cut` and `awk` to make the name of the second read. This then gets sent (first) to `TrimGalore` to trim reads. 

A note, because here we are doing this in markdown format (in R-studio), we can't use the workspace to keep variable names. We therefore need to redefine variables in each code chunk... This tends to be limited to file paths and the like. 

An important note: when I do this, I usually run the for-loop without the calls to `TrimGalore` in order to check that I have formatted `cut` and `awk` properly to generate the names that I like. I'll walk through this below.

```{bash, eval = FALSE}
basedir=/Volumes/BlytheLab_Files/HTSeq/Shelby/200824_cg_seChIP_pilot
rawdir=Raw_Data
trimdir=Trimmed_Reads
source /Applications/TrimGaloreEnv/bin/activate
for file in $basedir/$rawdir/*_R1_001.fastq.gz
do
filename=$(basename "$file")
echo "path to the first read"
echo $basedir/$rawdir/$filename
echo "base name for the first read"
echo $filename
filename2=$(echo $filename |cut -d_ -f1-3|awk '{print $1"_R2_001.fastq.gz"}')
echo "path to the second read"
echo $basedir/$rawdir/$filename2
echo "base name for the second read"
echo $filename2
echo
echo
# trim_galore -o $basedir/$trimdir --gzip --fastqc --paired $basedir/$rawdir/$filename $basedir/$rawdir/$filename2
done
deactivate
```

If everything looks good, then you can uncomment the line that runs `TrimGalore`. 

```{bash, eval = FALSE}
basedir=/Volumes/BlytheLab_Files/HTSeq/Shelby/200824_cg_seChIP_pilot
rawdir=Raw_Data
trimdir=Trimmed_Reads
source /Applications/TrimGaloreEnv/bin/activate
for file in $basedir/$rawdir/*_R1_001.fastq.gz
do
filename=$(basename "$file")
echo "path to the first read"
echo $basedir/$rawdir/$filename
echo "base name for the first read"
echo $filename
filename2=$(echo $filename |cut -d_ -f1-3|awk '{print $1"_R2_001.fastq.gz"}')
echo "path to the second read"
echo $basedir/$rawdir/$filename2
echo "base name for the second read"
echo $filename2
echo
echo
trim_galore -o $basedir/$trimdir --gzip --fastqc --paired $basedir/$rawdir/$filename $basedir/$rawdir/$filename2
done
deactivate
```

This should run and it will take a while to perform adapter trimming. After this is complete, we now need to perform mapping using bowtie2. Again, the exact details of the code will vary depending on whether this is a paired-end or single-end sequencing dataset. 

What the code below does is to perform mapping using near default parameters for `bowtie2`. It then pipes the output to `samtools view` (while also saving the mapping log (percent of reads mapped successfully, et cetera)). Next it uses `samtools sort` to sort the reads. Then it uses `picard CleanSam` to tidy up the file, and `picard MarkDuplicates` to mark, but not remove, suspected duplicates. Finally, it uses `samtools index` to create the necessary index file for the .bam file. 

Some notes about running this yourself: this code uses parallel processing, and I have this set to use up to 10 cores. Lab computers have only 4 available cores, and you need to change the following parameters to avoid an error:

in the call to `bowtie2`, change `-p 10` to `-p 4`. Also, in the call to `samtools sort`, change parameter `-@ 8` to `-@ 4`. If there are additional performance issues with `samtools sort`, the parameter `-m 2G` tells the program to allocate 2 gigabytes of memory to the sorting process, you may need to reduce this, but there should be sufficient memory on Lab computers to handle this.

Finally, the line at the top beginning `index=...` points `bowtie2` to precomputed indices for the genome. The path below is good for my office computer, but will need to be changed on other systems. This directory is located (currently) in Dropbox, but we may move it to the NAS if that makes more sense to do.

```{bash, eval = FALSE}
basedir=/Volumes/BlytheLab_Files/HTSeq/Shelby/200824_cg_seChIP_pilot
trimdir=Trimmed_Reads
mapdir=Mapped_Reads
index=/Users/sblythe/Dropbox/BlytheLab_HTSeq/Bowtie_Indices/dm6
for file in $basedir/$trimdir/*R1_001_val_1.fq.gz
do
filename=$(basename "$file")
echo "path to the first trimmed read"
echo $basedir$trimdir/$filename
echo "base name for the first trimmed read"
echo $filename
filename2=$(echo $filename |cut -d_ -f1-3|awk '{print $1"_R2_001_val_2.fq.gz"}')
echo "base name for the second trimmed read"
echo $filename2
shortname=$(echo $filename |cut -d_ -f1-2)
echo "short name for the output files"
echo $shortname
bowtie2 -p 10 -X 2000 -x $index -1 $basedir/$trimdir/$filename -2 $basedir/$trimdir/$filename2 2> $basedir/$mapdir/$shortname.bowtie.mapping.log.txt | \
samtools view -bS - > $basedir/$mapdir/$shortname.mapped.bam
samtools sort -@ 8 -m 2G -T $basedir/$mapdir/temp -o $basedir/$mapdir/temp1.bam $basedir/$mapdir/$shortname.mapped.bam
java -jar /Applications/picard.jar CleanSam \
I=$basedir/$mapdir/temp1.bam \
O=$basedir/$mapdir/temp2.bam
java -jar /Applications/picard.jar MarkDuplicates \
I=$basedir/$mapdir/temp2.bam \
O=$basedir/$mapdir/$shortname.mapped.md.bam \
M=$basedir/$mapdir/$shortname.dup.metrics.txt
rm $basedir/$mapdir/temp1.bam $basedir/$mapdir/temp2.bam $basedir/$mapdir/$shortname.mapped.bam
samtools index $basedir/$mapdir/$shortname.mapped.md.bam
done
