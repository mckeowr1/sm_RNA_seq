library(tidyverse)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)




files <- list.files()


test <- files[1]


df <- data.table::fread() 



#Plot the average number of reads per gene

#Pull Gene name and Average number of reads out of file names
df <- as.data.frame(files) %>%
  mutate(num_reads =  str_extract(.$files, pattern = "_[:digit:]+" )) %>%  #Pull out average number of reads per bin from file name
  mutate(num_reads = str_replace(.$num_reads, pattern = "_", "")) %>%  #Clean up underscore
  mutate(num_reads = as.numeric(num_reads)) %>% 
  mutate(num_reads_log10 = log(num_reads, 10)) %>% 
  mutate(gene_id = str_extract(.$files, pattern = "FBgn[:digit:]+"))

genes <- as.tibble(genes(txdb)) %>% #Get the genes as a gRanges obj
  dplyr::select(gene_id)

all_gene_df <- left_join(genes, df, by = "gene_id") %>% 
  mutate(num_reads = replace_na(num_reads, 0)) #Replace NA with 0

#Density Plot
ggplot(all_gene_df, aes(x=num_reads)) + 
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8)

#Histogram
ggplot(all_gene_df, aes(x=num_reads)) + 
  geom_histogram(fill="#69b3a2", color="#e9ecef", alpha=0.8) + 
  labs(x = "Average Number of Histone Cluster Reads Throughout Gene", 
       y = " Number of Genes"
      )


#Average Number of Reads for Multi_mapping genes 

ggplot(df, aes(x = num_reads)) + 
  geom_histogram(fill="#69b3a2", color="#e9ecef", alpha=0.8, bins = 100)

ggplot(df, aes(x = num_reads_log10)) + 
  geom_histogram(fill="#69b3a2", color="#e9ecef", alpha=0.8, bins = 100) + 
  labs(x = "Log10 Average Number of Histone Cluster Reads Throughout Gene", 
       y = " Number of Genes", 
       )


## Add some expression data ## 




## Annotate Genes with Flybase Biologic Process ##

#Write a file that can be input into flybase
multi_genes <- df$gene_id

rand_genes <- sample_n(genes, 1446) %>%  
  .$gene_id

writeLines(multi_genes, "multi_mapped_genes.txt")
writeLines(rand_genes, "rand_genes.txt")
#Load in the Flybase Biological Process data 

multi_mapped_bio <- data.table::fread("20211117_multimapped_genes_biological_process.txt") %>%  
  rename("gene_id" = "#SUBMITTED ID")

rand_bio <- data.table::fread("20211117_rand_genes_biological_process_2.txt") %>%  
  rename("gene_id" = "#SUBMITTED ID")

#Create Column for Biological Process of Interest 

multi_mapped_key_process <- multi_mapped_bio %>%
  mutate(unannotated = str_detect(.$GO_BIOLOGICAL_PROCESS, "^-")) %>% #Mark genes with no annotation  
  mutate(transcription = str_detect(.$GO_BIOLOGICAL_PROCESS, "transcription")) %>%  
  mutate(expression = str_detect(.$GO_BIOLOGICAL_PROCESS, "expression")) %>%  
  mutate(embryonic = str_detect(.$GO_BIOLOGICAL_PROCESS, "embryo")) %>%  
  mutate(oogenesis = str_detect(.$GO_BIOLOGICAL_PROCESS, "oogenesis"))
  

multi_gene_data <- left_join(df ,multi_mapped_key_process, by = c("gene_id")) %>%  
  dplyr::select("gene_id", "num_reads", "unannotated", "expression", "transcription", "embryonic", "oogenesis")



multi_summary_true <- multi_gene_data %>% 
  summarise_at(vars(unannotated:oogenesis), mean) %>%  #mean will tell us the number of true records
  gather() %>%  
  mutate(percent_multi_mapped_genes = value * 100 )

rand_key_process <- rand_bio %>%
  mutate(unannotated = str_detect(.$GO_BIOLOGICAL_PROCESS, "^-")) %>% #Mark genes with no annotation  
  mutate(transcription = str_detect(.$GO_BIOLOGICAL_PROCESS, "transcription")) %>%  
  mutate(expression = str_detect(.$GO_BIOLOGICAL_PROCESS, "expression")) %>%  
  mutate(embryonic = str_detect(.$GO_BIOLOGICAL_PROCESS, "embryo")) %>%  
  mutate(oogenesis = str_detect(.$GO_BIOLOGICAL_PROCESS, "oogenesis"))

rand_gene_data <- rand_key_process %>% 
  dplyr::select("gene_id", "unannotated", "expression", "transcription", "embryonic", "oogenesis")

rand_summary_true <- rand_gene_data %>%  
  summarise_at(vars(unannotated:oogenesis), mean) %>%  #mean will tell us the number of true records
  gather() %>%  
  mutate(percent_multi_mapped_genes = value * 100 )

#Plots for all Genes 
ggplot(multi_summary_true, aes(x = "", y=value, fill = key)) + 
  geom_bar(stat = "identity") + 
  coord_polar("y", start = 0) +
  theme_void() + 
  labs(title = "Biologic Process Annotation For Histone Cluster Multimapped Genes, N = 1446")

ggplot(rand_summary_true, aes(x = "", y=value, fill = key)) + 
  geom_bar(stat = "identity") + 
  coord_polar("y", start = 0)+ 
  theme_void() +
  labs(title = "Biologic Process Annotation For Randomly Selected Genes, N = 1446")





#Plots for the highest mapping genes
multi_gene_data %>% 
  dplyr::filter(num_reads > 10000) 


## Read Distribution along genes with multi_mappings ## 

#Blank Data frame to store the data from the loop
gene_percent_mapping <- data.frame()

for(f in files){ 
  
  df <- data.table::fread(f)

  gene_id <- str_extract(f, pattern = "FBgn[:digit:]+")

  #bin_size <- 100 # 100 bp bins set in gene_mapping.R

  num_bins <- max(df$Gene_Bin)


  bins_w_reads <- df %>% 
    dplyr::filter(num_reads > 0) %>%  
    group_by(Gene_Bin) %>%  
    summarise(number_reads = sum(num_reads)) %>% #If reads that map from different bins of the ROI map to the same region of the gene sum the num of reads
    dplyr::mutate(bin_percent = (Gene_Bin / num_bins) * 100 ) %>%
    dplyr::mutate(gene_id = gene_id) %>%
    dplyr::select(gene_id, number_reads, bin_percent)
  
  gene_percent_mapping <- bind_rows(gene_percent_mapping, bins_w_reads)
}


#Scatter Plot of number of reads by bin percent
ggplot(gene_percent_mapping, aes(x = bin_percent, y = number_reads)) + 
  geom_point() + 
  labs(x = "Bin Number/Total Bins per Gene", 
       y = "Number of Reads Per Bin Location", 
       title = "Number of Reads Per Gene Bin for Histone Cluster Multimapped Genes") 


#Log 10 scale - normal scale is fine
# gene_percent_mapping %>% 
#   dplyr::mutate(number_reads_log10 = log(number_reads, 10 )) %>% 
#   ggplot(aes(x = bin_percent, y = number_reads_log10)) + 
#   geom_point()








