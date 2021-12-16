library(tidyverse)

#This is a test script that I am hoping to develop into a way to filter through high mapping regions by looking at read patterns
#Takes all the high mapping genomic bins (over some threshold) and creates a bed file for the genomic loci 
#That bed file can then be used to subset the reads that map from one bin to a particular loci using 


setwd("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/analysis")

#Load Binned Genome DF
bins_df <- data.table::fread("binned_dm6_100000bins.tsv")

#Load Multi Mapping DF 
histone <- data.table::fread("histone_cluster_meltedmatrix.tsv")


#Filter for High mapping reads
high_mapping <- histone %>% 
  dplyr::filter(Log10_ReadCount > 5) %>%  
  dplyr::filter(GenomeBin != 215) %>% # Remove Anything that maps to histone cluster
  dplyr::filter(GenomeBin != 216)



#Make a list of bed files to loop through



for (row in 1:nrow(high_mapping)) {
  df <- high_mapping[row] 
  bed <- as.character(df[,2])
  genomic_bin <- as.numeric(df[,1])
  
  file_name <- str_c(genomic_bin, bed, sep = "_")
  
  ##Write out a bed file##
  df_bin <- bins_df %>%
    dplyr::filter(bin_id == genomic_bin) %>%  
    select(seqnames, start, end)
  
  
  setwd("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/histone_cluster/high_mappings/genomic_beds/")
  
  #Write out a file with range of the genomic bed of interest
  write_tsv(df_bin, glue::glue("{file_name}.bed"), col_names = FALSE)
  
  
  print(df_bin)
}


test_bed <- bins_df %>% dplyr::filter(bin_id == 65) %>%  
  select(seqnames, start, end)


write_tsv(test_bed, "test_region.bed", col_names = FALSE)

test_bed 




# #Check Mapping to histone cluster 
# 
# histone_mappers <- histone %>%  
#   dplyr::filter(GenomeBin == 215) %>%  
#   dplyr::filter(Log10_ReadCount > 0)
# 
# two_sixty <- histone %>%
#   dplyr::filter(GenomeBin == 260) %>%  
#   dplyr::filter(Log10_ReadCount > 0)



