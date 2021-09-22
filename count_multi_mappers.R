library(data.table)
library(tidyverse)

df <- data.table::fread("SRR1187947_trimmed.fq.gz.bowtie.multimapping_non_primary.sam") %>%  
  group_by(V1) %>% 
  summarise(n = n(), 
            min_qual = min(V5), 
            max_qual = max(V5))
hist(df$n, xlim = c(10000
                    ,16000))


