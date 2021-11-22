library(GenomicAlignments)
library(bedr)


#Read-In Bed with bedr

#603 Lines 0.06
system.time(bed_to_granges("bin10_bed_ac"))

normal_file <- data.table::fread("bin10_bed_ac")

large_file <- rbind(normal_file, normal_file, normal_file, normal_file)

write_tsv(large_file, "large_test_bed")

fast_read <- function(file){ 

df <- data.table::fread(file)


gr <- makeGRangesFromDataFrame(df, 
                         seqnames.field = "V1", 
                         start.field = "V2", 
                         end.field = "V3", 
                         
                         
                         )
}

#Faster 0.016 seconds ! 603 enteries 
system.time(fast_read("bin10_bed_ac"))


#Time for a larger file (2412 Obs)
system.time(fast_read("large_test_bed"))

# 
# granges(
#   seqnames = df$V1
#   ranges = I
# )
