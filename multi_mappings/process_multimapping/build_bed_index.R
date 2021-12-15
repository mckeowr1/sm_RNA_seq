library(data.table)
library(tidyverse)

df <- data.table::fread("/path/to/bed_file") #Add path to bed file

build_bed_index <- function(bed, name){ 
  df <- rownames_to_column(bed) 

  read_key <- group_by(df, V4) %>%  
  summarise(index = paste0(as.character(rowname), collapse = ","))

  write_tsv(read_key, glue::glue("{name}.bedindex.tsv"))
}

build_bed_index(bed = df, name = ROI) 