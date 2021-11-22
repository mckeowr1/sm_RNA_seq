library(data.table)
library(tidyverse)

setwd("/Users/ryan/Documents/GitHub/sm_RNA_seq/sat_dna/")
read_count <- data.table::fread("count.out")

repeats <- unique(read_count$`#repeat`)

gaga_repeats<- str_subset(string = repeats, pattern = "GAGA")

gaga_reads <- read_count %>% dplyr::filter(`#repeat` %in% gaga_repeats)

gaga_elements_w_reads <- gaga_reads %>% dplyr::filter(`#count` > 0) %>% 
  .$`#repeat`

setwd("/Users/ryan/Documents/GitHub/sm_RNA_seq/sat_dna/ping_pong/")

writeLines(gaga_elements_w_reads, "gaga_elements_w_reads.txt")
