library(data.table)
library(tidyverse)

setwd("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/")


index <- data.table::fread("files/SRR1187947_mapped_verysensitive_local.mapped.bedindex.tsv")



#Filter the index based on Read Names Input

#Region

read_names <- readLines("files/small_readnames.txt")

for(read in read_names){ 
  get_bed_location(index, read)
  
  } 


#Function to Search through a Bed Index and Pull out line loactions 
get_bed_location <- function(bedindex, 
                             readname)
{
  ## At some point want to figure out how to run this for multiple reads## 
  # #Read in read names to character vector 
  # read_names <- readLines(read_name )
  
  read <- bedindex %>% dplyr::filter(V4 == readname) 
  
  lines <- as.character(read[,2]) %>%
    str_split(pattern = ",") %>%  #For some reason this converts to a list
    .[[1]] %>% #Make it not a list 
    sort()

    writeLines(lines, glue::glue("{readname}_lines.txt"))
  
  
} 



