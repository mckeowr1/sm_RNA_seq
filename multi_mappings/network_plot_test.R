library(data.table)
library(tidyverse)
library(igraph)
library(feather)
library(networkD3)


setwd("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/search_beds/")

index <- data.table::fread("SRR1187947_mapped_verysensitive_local.mapped_sorted.bedindex.tsv")

sm_index <- sample_n(index, 15)

r_index <- rename(index, read_name = V4)


#Read In the Big Bed file

#bed <- data.table::fread("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/search_beds/SRR1187947_mapped_verysensitive_local.mapped_sorted.bed", 
                  #select = c("V1", "V2", "V3"), nrows = 1000)


bedaa <-  data.table::fread("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/search_beds/SRR1187947_mapped_verysensitive_local_sortedaa.bed", 
                            select = c("V1", "V2", "V3")) %>%  
                      rownames_to_column() 
bedab <- data.table::fread("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/search_beds/SRR1187947_mapped_verysensitive_local_sortedab.bed", 
                           select = c("V1", "V2", "V3")) %>%  
  rownames_to_column() 

mu_aa <- bedaa %>% mutate(pos = str_c(V1, V2, V3, sep = "_"))
  
mu_bed_aa <- mu %>% select(rowname, pos) %>%  
  rename(index = rowname) %>%  
  mutate(index = as.numeric(index)) 


mu_ab <- bedab %>% 
  mutate(pos = str_c(V1, V2, V3, sep = "_")) %>% 
           select(rowname, pos) %>% 
           rename(index = rowname) %>%
           mutate(index = as.numeric(index))

write_feather(mu_bed, "mu_bedaa.feather")
write_feather(mu_ab, "mu_bedab.feather")

#Import a list of Readnames

setwd("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/42AB_test/binned_readnames")

reads <- readLines(small_readnames[1]) %>% as.data.table() %>%  
  rename(.,read_name = .)

#sm_reads <- sample_n(reads, 500) 

#Join the readnames to the index
setkey(reads, .)
setkey(index, V4)


new_test_df <- merge(reads, r_index, by = "read_name", all.x = TRUE) 


df3 <- new_test_df %>% 
  separate_rows(index, sep = ",") %>%  
  mutate(index = as.numeric(index)) %>% 
  dplyr::filter(between(index, 100000001, 200000000)) %>%  
  mutate(index = index - 100000000)

156682524
100000000


#Merge the Reads + Index DF to the Bed File 

merged <- merge(df3, mu_ab, by = "index", all.x = TRUE)

plot_df <- merged %>%  
  select(read_name, pos)

set.seed(1492)

l <- layout.fruchterman.reingold(g, niter=5000, area=vcount(ig)^4*10)


g <- graph_from_data_frame(plot_df, directed = FALSE)
l <- layout.fruchterman.reingold(g, niter=5000, area=vcount(g)^4*10)
l2 <- layout_nicely(g, dim=2)
plot.igraph(g, 
            layout = l2, 
            edge.arrow.size=0.5, 
            vertex.label=NA,
            vertex.shape="circle", 
            vertex.size=1, 
            vertex.label.color="black", 
            edge.width=0.5)

p <- simpleNetwork(plot_df, height ="500px", width = "500px", 
                   Source = 1,
                   Target = 2, 
                   linkDistance = 10, 
                   charge = -900, 
                   fontSize = 1, 
                   zoom = T
                   )
p


saveWidget(p, "1bin_interactivenetwork.html")



