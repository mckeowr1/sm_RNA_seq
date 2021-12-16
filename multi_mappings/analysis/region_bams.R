library(data.table)
library(tidyverse)

# setwd("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/analysis/fbgn_0027339/")

# allreads = readGAlignments(file, 
#                            param = ScanBamParam(flag = scanBamFlag(isSecondaryAlignment = FALSE, isUnmappedQuery = FALSE, isDuplicate = FALSE), 
#                                                 what = c('flag','mrnm','mpos','mapq','isize')),
#                            use.names = FALSE)

# gene_reads <- readGAlignments("fbgn_0027339.sam", use.names = TRUE) 

# check <- as.data.frame(gene_reads)

bin_readnames <- data.table::fread("fbgn_0027339_histone_reads.bed", select= c("V4"))%>% 
  .$V4



writeLines(bin_readnames, "fbgn_0027339_histone_readnames.txt")


check <- readGAlignments("fbn_0027339_histone.bam")
