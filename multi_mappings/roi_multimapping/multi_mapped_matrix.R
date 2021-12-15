library(GenomicAlignments)
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(data.table)






multi_map_matrix <- function(roi_name, genome_bin_size = 100000, bed_dir, out_dir){
    #Load the reference genome 
    bins = tileGenome(seqlengths(Dmelanogaster), tilewidth = genome_bin_size, cut.last.tile.in.chrom = TRUE)
    bins = bins[seqnames(bins) %in% paste0('chr', c('2L','2R','3L','3R','4','X'))] # this can be commented out if you want to keep all chromosomes/contigs.
    bins = bins[width(bins) == genome_bin_size] 

    #Set WD to the bed files for the roi
    setwd(bed_dir)
    
    #Get a list of the bins 
    bin_beds <- list.files()

    #Order the beds for naming purposes 
    library(stringr)
    ordered_beds = stringr::str_sort(bin_beds, numeric = T)

    
    #Define function to make a vector from a bed files intersection with the genome 

    make_vector<- function(file){ 
        print(file)
  
         df <- data.table::fread(file) #Read in the bedfile
  
  
        gr <- makeGRangesFromDataFrame(df, 
                                 seqnames.field = "V1", 
                                 start.field = "V2", 
                                 end.field = "V3", 
                                 strand.field = "V6")
         c <- countOverlaps(bins, gr)
  
        return(c)                             

                    }
    
    
    #Apply the Vector Function on all the Bin Beds
    myfiles = lapply(ordered_beds, make_vector)

matrix = do.call("cbind", myfiles)

# Melt the matrix
library(reshape2)
library(tidyverse)

longdata <- reshape2::melt(matrix) %>%
    rename(GenomeBin = Var1, HistoneClusterBin = Var2, ReadCount = value) %>%
    mutate(Log10_ReadCount = log(value + 1, 10))


write_tsv(longdata, glue::glue("{out_dir}/{roi_name}_multimapped_matrix.tsv")


}





#Load the mack