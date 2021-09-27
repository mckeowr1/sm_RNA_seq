library(data.table)
library(tidyverse)

setwd("/Users/ryan/Documents/GitHub/sm_RNA_seq/SRR1187947/Multi_Mapped_Reads/")

df <- data.table::fread("SRR1187947_trimmed.fq.gz.bowtie.multimapping_non_primary.sam") #%>%  
  #group_by(V1) %>% 
  #summarise(n = n(), 
            #min_qual = min(V5), 
            #max_qual = max(V5))


df %>%  
  dplyr::filter(V12 != V13)



m1 <- mean(df$n)
median(df$n)
max(df$n)

glue::glue("The mean # of possible mappings is {m1}")

#Qualtiy Score check


df %>%  
ggplot(aes(x = n)) + 
  geom_freqpoly()




