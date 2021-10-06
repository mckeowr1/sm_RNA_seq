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



