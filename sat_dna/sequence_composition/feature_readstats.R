library(seqinr)
library(ggseqlogo)
library(ggplot2)
library(Biostrings)
library(seqTools)

gagaca<- read.fasta("/Users/ryan/Documents/GitHub/sm_RNA_seq/sat_dna/GAGACA_trimmed.fasta")

length(gagaca)

tenth_base <- function(x){
  return(x[10])
}

table(mapply(tenth_base, gagaca))


#write a function to return the first base
first_base <- function(x){
  return(x[1])
}

table(mapply(first_base, gagaca))


lengths <- table(getLength(gagaca))
barplot(lengths)


readDNAStringSet("/Users/ryan/Documents/GitHub/sm_RNA_seq/sat_dna/GAGACA_trimmed.fasta", "fasta")
# afmc=consensusMatrix(f1, baseOnly=T,as.prob = T)
# tafmc = t(afmc)

library(seqTools)
# Reads fastq file
fq=fastqq("/Users/ryan/Documents/GitHub/sm_RNA_seq/sat_dna/GAGACA_trimmed.fq")
# Plots nucleotide frequency
plotNucFreq(fq,1)


rsp <- read.fasta("/Users/ryan/Documents/GitHub/sm_RNA_seq/sat_dna/Rsp_trimmed.fasta")

#piRNAs usually have a U at the first base 
table(mapply(first_base, rsp))
#piRNAs usually have a A at the 10th base 
table(mapply(tenth_base, rsp))




lengths <- table(getLength(rsp))
barplot(lengths)

#ggplot()+


