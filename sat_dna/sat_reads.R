library(GenomicAlignments)
library(bedr)



setwd("/Users/ryan/Documents/GitHub/sm_RNA_seq/sat_dna/")


allreads = readGAlignments("SRR1187947_mapped_verysensitive_local_wei2020_bestreads.mapped.sorted.bam", 
                           param = ScanBamParam(flag = scanBamFlag(isSecondaryAlignment = FALSE, isUnmappedQuery = FALSE, isDuplicate = FALSE), 
                                                what = c('flag','mrnm','mpos','mapq','isize')),
                           use.names = FALSE)



#Filter for reads that map to (AAGAG)n




#gaga_rna_sats <- bed_granges("gaga_sats/gaga_rna_sites_karpen2019.bed")

gff_granges <- function(file){
  t_df <- data.table::fread(file)
  
  makeGRangesFromDataFrame(t_df, 
                           seqnames.field = "V1",
                           start.field = "V4", 
                           end.field = "V5", 
                           strand.field = "V7", 
                           keep.extra.columns = TRUE)
}




aagag_sats <- gff_granges("/Users/ryan/Documents/GitHub/sm_RNA_seq/sat_dna/aagag_repeats.gff3")


aagaga_reads <- subsetByOverlaps(allreads, aagag_sats)


test_df <- as.data.table(aagaga_reads)

hist(test_df$width)







as.data.frame(aagaga_reads)

