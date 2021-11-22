library(tidyverse)

#setwd("/Users/ryan/Documents/GitHub/sm_RNA_seq/sat_dna/ping_pong/")

#features  <- readLines("gaga_elements_w_reads.txt")

setwd("/Users/ryan/Documents/GitHub/sm_RNA_seq/sat_dna/ping_pong/rsp_piProfiles/")

features <- list.files()

features[1]

df <- data.table::fread(features[1], col.names = c("overlap", "n_pairs", "z-score", "probability"))


#Check the n_pairs for any actually overlapping reads 

check <-  df %>% filter(n_pairs > 1)

features_w_overlapping_reads <- c()
for(f in features){ 

  element_name <-stringr::str_extract(f, ".+(?=[:punct:])")

  df <- data.table::fread(f, col.names = c("overlap", "n_pairs", "z-score", "probability")) %>%  
    filter(n_pairs > 1)
  
if(nrow(df) == 0 ){ 
  cat(glue::glue("There are no overlapping reads for {element_name}"), "\n" )
  } else 
    features_w_overlapping_reads<- append(features_w_overlapping_reads, f)

} 


#Change Element Name!!
setwd("/Users/ryan/Documents/GitHub/sm_RNA_seq/sat_dna/ping_pong/sat_piProfiles/")
pdf("bp_sat.pdf")
for(f in features_w_overlapping_reads){ 
  element_name <-stringr::str_extract(f, ".+(?=[:punct:])")
  
   df <- data.table::fread(f, col.names = c("overlap", "n_pairs", "z-score", "probability"))
  
  f <- df %>% 
    ggplot(aes(x = overlap, y = n_pairs)) + 
    geom_bar(stat = "identity") + 
    ggtitle(glue::glue("{element_name}"))
  
  plot(f)

} 
dev.off()

plot(f)

  
#stringr::str_extract(features[1], ".+(?=[:punct:])")
#print(features[1])

