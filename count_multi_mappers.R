library(data.table)
library(tidyverse)

df <- data.table::fread("multi_mapped.sam") %>%  
  group_by(V1) %>% 
  summarise(n = n(), 
            min_qual = min(V5), 
            max_qual = max(V5))
