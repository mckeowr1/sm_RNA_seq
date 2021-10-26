library(GenomicAlignments)
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(data.table)

bins = tileGenome(seqlengths(Dmelanogaster), tilewidth = 100000, cut.last.tile.in.chrom = TRUE)
bins = bins[seqnames(bins) %in% paste0('chr', c('2L','2R','3L','3R','4','X'))] # this can be commented out if you want to keep all chromosomes/contigs.
bins = bins[width(bins) == 100000] 

make_vector<- function(file){ 
  print(file)
  
  df <- data.table::fread(file)
  
  
  gr <- makeGRangesFromDataFrame(df, 
                                 seqnames.field = "V1", 
                                 start.field = "V2", 
                                 end.field = "V3")
  c <- countOverlaps(bins, gr)
  
  return(c)                             

}


setwd("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/42AB_test/42AB_beds/")
bin_beds <- list.files()

library(stringr)
ordered_beds = str_sort(bin_beds, numeric = T)


myfiles = lapply(ordered_beds, make_vector)

names(myfiles) <- ordered_beds

matrix = do.call("cbind", myfiles)



heatmap(matrix, Colv = NA, Rowv = NA, na.rm = FALSE)



#GGplot a Matrix

library(reshape2)
library(ggplot2)

longdata <- melt(matrix)

ggplot(longdata, aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() + 
  scale_fill_gradient(low= "grey90", high="red")

#Overlap the files with Binned References


#Annotation Data

forty=GRanges(seqnames =c("chr2R"), 
              ranges=IRanges(start =c(6255432), end=c(6499291)))
forty_matrix = countOverlaps(bins, forty)


setwd("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/annotation_data/")

library(tidyverse)
ext_high_ov <- data.table::fread("extremly_highovary.txt") %>%  
  mutate(LOCATION_ARM = str_c("chr", .$LOCATION_ARM)) #Add chr to the location arm so that it can be intersected
  
hi_ov_gr <- makeGRangesFromDataFrame(high_ov, 
                         seqnames.field = "LOCATION_ARM", 
                         start.field = "LOCATION_MIN", 
                         end.field = "LOCATION_MAX")



v_high_ov <- data.table::fread("Veryhigh_ovary_expression.tsv") %>% 
  mutate(LOCATION_ARM = str_c("chr", .$LOCATION_ARM))


v_hi_ov_gr <- makeGRangesFromDataFrame(high_ov, 
                                       seqnames.field = "LOCATION_ARM", 
                                       start.field = "LOCATION_MIN", 
                                       end.field = "LOCATION_MAX")


get_expression_vec <- function(file){
  df <- data.table::fread(file) %>%  
    mutate(LOCATION_ARM = str_c("chr", .$LOCATION_ARM))
  
  gr <- makeGRangesFromDataFrame(df, 
                                 seqnames.field = "LOCATION_ARM", 
                                 start.field = "LOCATION_MIN", 
                                 end.field = "LOCATION_MAX")
  
  bin_IDs<- countOverlaps(bins, gr) %>%  
    as.data.frame() %>% 
    rownames_to_column %>%  
    filter(. == 1) %>% 
    pull(rowname)
  
  return(bin_IDs)
}


vhi_ovary <- get_expression_vec("Veryhigh_ovary_expression.tsv")
hi_ovary <- get_expression_vec("high_ovary.tsv")
vlow_ovary <- get_expression_vec("verylow_ovary.tsv")

#Make a character vector of bins with 42AB
forty_bins<- forty_matrix %>% as.data.frame() %>% rownames_to_column() %>% filter(. == 1) %>%  
  pull(rowname)

mutate_cond <- function(.data, condition, ..., envir = parent.frame()) {
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data[condition, ] %>% mutate(...)
  .data
}



df<- countOverlaps(bins, hi_ov_gr) %>%  as.data.frame() %>%  
  rownames_to_column() %>%  
  mutate(color = ifelse( . >= 1, "#6a040f", NA)) %>%  
  mutate_cond(rowname %in% forty_bins, color = "black") %>%  
  #mutate_cond(rowname %in% vhi_ovary, color = "#9d0208") %>%  
  #mutate_cond(rowname %in% high_ovary, color = "#d00000" ) 
  mutate_cond(rowname %in% vlow_ovary, color = "#abc4ff")
  

colors <- df$color  


setwd("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/42AB_test/plots")

pdf("42AB_multimapping_profile.pdf")

heatmap(matrix, Colv = NA, Rowv = NA, RowSideColors = colors , xlab = "42AB 100 bp Bin", ylab = "Genomic Pos 100000 bp Bin")
dev.off()


