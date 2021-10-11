library(GenomicAlignments)
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(data.table)

bins = tileGenome(seqlengths(Dmelanogaster), tilewidth = 100000, cut.last.tile.in.chrom = TRUE)
bins = bins[seqnames(bins) %in% paste0('chr', c('2L','2R','3L','3R','4','X'))] # this can be commented out if you want to keep all chromosomes/contigs.
bins = bins[width(bins) == 100000] 

make_vector<- function(file){ 
  
  df <- data.table::fread(file)
  
  
  gr <- makeGRangesFromDataFrame(df, 
                                 seqnames.field = "V1", 
                                 start.field = "V2", 
                                 end.field = "V3")
  c <- countOverlaps(bins, gr)
  
  return(c)                             

}


setwd("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/42AB_test/bin_beds")
bin_beds <- list.files()

library(stringr)
ordered_beds = str_sort(bin_beds, numeric = T)


myfiles = lapply(ordered_beds, make_vector)



matrix = do.call("cbind", myfiles)

heatmap(matrix, Colv = NA, Rowv = NA)

#Overlap the files with Binned References


#Annotation Data

forty=GRanges(seqnames =c("chr2R"), 
              ranges=IRanges(start =c(6255432), end=c(6499291)))
forty_matrix = countOverlaps(bins, forty)


setwd("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/annotation_data/")

library(tidyverse)
high_ov <- data.table::fread("extremly_highovary.txt") %>%  
  mutate(LOCATION_ARM = str_c("chr", high_ov$LOCATION_ARM)) #Add chr to the location arm so that it can be intersected
  

hi_ov_gr <- makeGRangesFromDataFrame(high_ov, 
                         seqnames.field = "LOCATION_ARM", 
                         start.field = "LOCATION_MIN", 
                         end.field = "LOCATION_MAX")

df <- countOverlaps(bins, hi_ov_gr) %>%  as.data.frame() %>%  
  rownames_to_column() %>%  
  mutate(color = ifelse( . >= 1, "#FF0000", "#FFFF00")) 

colors <- df$color  


heatmap(matrix, Colv = NA, Rowv = NA, RowSideColors = colors)



