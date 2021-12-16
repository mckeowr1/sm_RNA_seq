


### Get Histone Cluster H4 Gene Bin ID ##


## Set up the Reference Genome ## 
bins = tileGenome(seqlengths(Dmelanogaster), tilewidth = 100000, cut.last.tile.in.chrom = TRUE)
bins = bins[seqnames(bins) %in% paste0('chr', c('2L','2R','3L','3R','4','X'))] # this can be commented out if you want to keep all chromosomes/contigs.
bins = bins[width(bins) == 100000] 

roi=GRanges(seqnames =c("chr2L"), 
            ranges=IRanges(start =c(21415940), end=c(21543673))) 

## Subset the Region of Interest from the tiled Genome 

roi_t <- subsetByOverlaps(bins, roi)

#Give the Bins an id 
bins$binID = 1:length(bins)


gene_of_interest=GRanges(seqnames =c("chr2R"), 
                         ranges=IRanges(start =c(10977956), end=c(11116378))) 


subsetByOverlaps(bins, gene_of_interest)


21,529,242 21,529,769


#21,534,085 21,534,610

## Get the Readnames list to filter down the Bam file 


setwd("~/Documents/GitHub/sm_RNA_seq/multi_mappings/analysis/H4")

readnames <- data.table::fread("h4_histone_reads.bed", select = c("V4")) %>% 
  .$V4

writeLines(readnames, "h4_histone_reads.txt")





