library(data.table)

setwd("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/search_beds/")


all_bed <- data.table::fread("SRR1187947_mapped_verysensitive_local.mapped_sorted.bed")


aa <- data.table::fread("SRR1187947_mapped_verysensitive_local_sortedaa.bed", select = c("V1", "V2", "V3"), nThread = 6)
ab <- data.table::fread("SRR1187947_mapped_verysensitive_local_sortedab.bed", select = c("V1", "V2", "V3"), nThread = 6 )
ac <- data.table::fread("SRR1187947_mapped_verysensitive_local_sortedac.bed", select = c("V1", "V2", "V3"), nThread = 6)
ad <- data.table::fread("SRR1187947_mapped_verysensitive_local_sortedad.bed", select = c("V1", "V2", "V3"), nThread = 6)
ae <- data.table::fread("SRR1187947_mapped_verysensitive_local_sortedae.bed", select = c("V1", "V2", "V3"), nThread = 6)
af <- data.table::fread("SRR1187947_mapped_verysensitive_local_sortedaf.bed", select = c("V1", "V2", "V3"), nThread = 6)
ag <- data.table::fread("SRR1187947_mapped_verysensitive_local_sortedag.bed", select = c("V1", "V2", "V3"), nThread = 6)
ah <- data.table::fread("SRR1187947_mapped_verysensitive_local_sortedah.bed", select = c("V1", "V2", "V3"), nThread = 6)






sm_index <- sample_n(index, 100 )

#List only search beds
search_beds = list.files(pattern = "a..bed")





