library(data.table)
library(tidyverse)




df <- data.table::fread("SRR1187947_mapped_verysensitive_local.mapped.bed", select = c("V4"))  

df2 <- rownames_to_column(df)


df3 <- group_by(df2, V4) %>%  
  summarise(index = paste0(as.character(rowname), collapse = ","))


dplyr::filter(df3, V4 == "SRR1187947.9999995") 
#Each read that is aligned is accounted for! 
write_tsv(df3, "SRR1187947_mapped_verysensitive_local.mapped.bedindex.tsv")
