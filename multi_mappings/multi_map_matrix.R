library(GenomicAlignments)
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(bedr)
dm6 = Dmelanogaster


setwd("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/42AB")

#### Setup the Reference Genome & The ROI Genome ####

bins = tileGenome(seqlengths(Dmelanogaster), tilewidth = 100, cut.last.tile.in.chrom = TRUE)
bins = bins[seqnames(bins) %in% paste0('chr', c('2L','2R','3L','3R','4','X'))] # this can be commented out if you want to keep all chromosomes/contigs.
bins = bins[width(bins) == 100] 


# Subset the Region of Interest From Tiled Genome # 

forty=GRanges(seqnames =c("chr2R"), 
              ranges=IRanges(start =c(6255432), end=c(6499291)))

#Get the tiled region of the genome
roi_t <- subsetByOverlaps(bins, forty)

#Name the bins 
roi_t$binID = 1:length(roi_t)

#### Assign Reads to BIN ####

#Read in Region BED 
roi_bed <- bed_to_granges("roi_bed/42AB.bed")

#sm_bed <- head(roi_bed, n= 1000)



## Deal with reads mapped to Region >1 ##

#Set the names to be Seq ID
names(roi_bed) <- mcols(roi_bed)[,"id"]

#Remove Duplicate Region Mappings
roi_bed_clean = roi_bed[unique(mcols(roi_bed)[,1])]


## Intersect the ROI_read_bed with the Tiled Genome Get a tiled ROI## 
overlaps <- mergeByOverlaps(roi_t, roi_bed_clean)

#### Get a Read Names File for each Bin #### 

library(tidyverse)


#
df <- overlaps[c(2, 4)] %>% as.data.frame()

reads_in_bin <- df %>% group_by(binID) %>% 
  summarize(reads = paste0(as.character(id), collapse = ",")) 


#Write this out to a file (name = bin_ID, Text = Read names)

setwd("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/42AB_test/binned_readnames/")
write_readnames_binned <- function(reads_in_bin) { 
  bin_name <- reads_in_bin[1]
  reads <- as.character(reads_in_bin[,2]) %>%  
    str_split(pattern = ",") %>% 
    .[[1]] %>%
    sort()
  
  writeLines(reads, glue::glue("bin{bin_name}_readnames.txt"))
  
} 

for(read in 1:nrow(reads_in_bin)){ 
  # print(reads_in_bin[read, 1:2])
  df <- reads_in_bin[read,]
  
  write_readnames_binned(df)

  
  }


test_file <- reads_in_bin[10,]
write_readnames_binned(test_file)




file_name <- test_file[1]
reads <- as.character(test_file[,2]) %>%  
  str_split(pattern = ",") 

writeLines(reads, glue::glue("bin{file_name}.txt"))

#write_tsv(reads_in_bin, "SRR1187947_mapped_verysensitive_local.mapped.bedindex.tsv")


writeLines(reads, "test")




