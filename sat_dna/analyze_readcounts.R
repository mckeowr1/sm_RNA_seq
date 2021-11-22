library(data.table)
library(tidyverse)

setwd("/Users/ryan/Documents/GitHub/sm_RNA_seq/sat_dna/")
read_count <- data.table::fread("count.out")

repeats <- unique(read_count$`#repeat`)


gaga_repeats<- str_subset(string = repeats, pattern = "GAGA")

gaga_reads <- read_count %>% dplyr::filter(`#repeat` %in% gaga_repeats)



ggplot(gaga_reads, aes( x= `#repeat`, y = `#count`))+ 
  geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 5))
  
